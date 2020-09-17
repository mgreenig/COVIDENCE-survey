import pandas as pd
import numpy as np
from scipy.stats import zscore, norm

from Map_survey_answers import q142_cleaned, q143_dosages, q143_units, survey_data
from Annotate_patients import drug_dictionary, bnf_classes, drug_classes, specific_drugs

# dictionary for patient drug classes
drug_class_doses = {}

# dictionary for patient drugs
specific_drug_doses = {}

# get DB ids for each survey answer
q142_db_ids = q142_cleaned.apply(lambda answer: drug_dictionary.get(answer) if drug_dictionary.get(answer) else set())

# function for aligning question indices in a survey answer mask to another survey answer series
def align_mask(mask, series):
    # get the question patterns from of the indices both the mask and answer series
    q1_pattern = mask.index[0][1][:5]
    q2_pattern = series.index[0][1][:5]
    # rename the question index in the mask
    mask_renamed = mask.rename({q_id: q_id.replace(q1_pattern, q2_pattern) for q_id in mask.index.get_level_values(1).unique()}, level = 1)
    # align it to the series
    _, mask_aligned = series.align(mask_renamed, fill_value = False)
    return mask_aligned

# function for standardising drug q143_dosages for a drugbank ID
def get_normalised_dosages(id):

    # if a single id is provided, convert into a set
    ids = set([id]) if isinstance(id, str) else id

    # make copy of the dosage question so the global object is not modified
    dosages = q143_dosages.copy()
    # filter for answers with the same DB id(s)
    drug_mask = q142_db_ids.apply(lambda drug_ids: True if ids == drug_ids else False)
    # change index name
    drug_mask = align_mask(drug_mask, q143_units)

    # if any drugs are hit, get the dosages
    if drug_mask.any():

        # get all dosage q143_units specified for the drug
        dosage_units = q143_units[drug_mask]
        mg_mask = align_mask(dosage_units == 1, dosages)
        micg_mask = align_mask(dosage_units == 2, dosages)
        invalid_unit_mask = align_mask((dosage_units != 1) & (dosage_units != 2), dosages)

        # get iqr and quantiles of the mg drug q143_dosages
        q1, q3 = dosages[mg_mask].quantile([0.25, 0.75])
        dose_iqr = q3 - q1

        # change the microgram values in the actual dosage values if the resulting answers are not outliers
        dosages[micg_mask] = np.where((dosages[micg_mask] / 1000 < q1-dose_iqr) | (dosages[micg_mask] / 1000 > q3+dose_iqr),
                                      dosages[micg_mask], dosages[micg_mask] / 1000)

        # get all dosage values with valid dosage units
        valid_unit_mask = mg_mask | micg_mask
        valid_dosages = dosages[valid_unit_mask]

        # calculate z score
        valid_dosages_scaled = zscore(valid_dosages)

        # normalised using probit function
        valid_dosages_norm = pd.Series(norm.cdf(valid_dosages_scaled), index = valid_dosages.index)

        # get invalid dosages
        invalid_dosages = dosages[invalid_unit_mask].apply(lambda dosage: 'NA')

        # combine the two dosage lists
        all_dosages = valid_dosages_norm.append(invalid_dosages)
        all_dosages = all_dosages.sort_index()

        # align to the larger series
        _, drug_dosages_aligned = dosages.align(all_dosages, fill_value=0)

    # otherwise save the dosages as a series of 0
    else:
        drug_dosages_aligned = pd.Series(0, index = dosages.index)

    return drug_dosages_aligned

# function for checking if a series of dosage answers consists of only NA and 0 values
is_na = lambda row: False if all(row == 0) else False if all(row != 'NA') else True

# function for combining dosage data across multiple drugs into a single value (either a sum or NA)
combine_func = lambda row: 'NA' if is_na(row) else row[row != 'NA'].sum()

# function for getting the dosage values for each class
def get_class_doses(drug_class):

    # get bnf entries that fall into the desired drug class
    drug_class_mask = bnf_classes['primary'].str.contains(drug_class) | bnf_classes['secondary'].str.contains(drug_class)
    class_drugs = bnf_classes[drug_class_mask]

    # isolate entries that in the drug class that correspond to single drugbank IDs, to avoid ambiguity with mixture products
    single_id_mask = class_drugs['db_id'].apply(lambda ids: len(ids) == 1)
    # get the db ids corresponding to the drug class
    class_db_ids = set().union(*class_drugs['db_id'][single_id_mask])

    # loop through ids in the class and add dosages to the dictionary of dosage data
    all_dosages = {}
    for id in class_db_ids:
        dosages = get_normalised_dosages(id)
        all_dosages[id] = dosages

    # make data frame from the dosage data for all drugs in the class
    dosage_df = pd.DataFrame(all_dosages)

    # group by patient, summing to get the total dose per patient within each class
    question_dosages = dosage_df.apply(combine_func, axis = 1)
    patient_dosages = question_dosages.groupby(level = 0).apply(combine_func)

    # align to the total set of survey answers
    _, patient_dosages_aligned = survey_data.align(patient_dosages, axis = 0, fill_value = 0)

    return patient_dosages_aligned

# function for getting the class
def get_drug_doses(drug):

    # get the drugbank id corresponding to the drug of interest
    drug_db_ids = drug_dictionary[drug]
    # get the normalised dosages
    dosages = get_normalised_dosages(drug_db_ids)
    # group by patient, summing to get the total dose per patient within each class
    patient_dosages = dosages.groupby(level = 0).apply(combine_func)
    # align to the total set of survey answers
    _, patient_dosages_aligned = survey_data.align(patient_dosages, axis = 0, fill_value = 0)

    return patient_dosages_aligned


if __name__ == '__main__':

    # label patient drug classes
    for drug_class in drug_classes:
        drug_class_doses[drug_class] = get_class_doses(drug_class)

    # label patient drugs
    for drug in specific_drugs:
        specific_drug_doses[drug] = get_drug_doses(drug)

    # put all the patient features into a single dictionary
    patient_dose_feature_dict = {**drug_class_doses, **specific_drug_doses}

    # make a data frame
    patient_dose_feature_df = pd.DataFrame(patient_dose_feature_dict)
    patient_dose_feature_df.insert(0, 'uid', survey_data['uid'])

    patient_dose_feature_df.rename({'^calcium$': 'calcium', 'oestrogens|androgens': 'sex hormone therapy',
                                    'antimuscarinics, other': 'antimuscarinics',
                                    'non-steroidal anti-inflammatory drugs': 'nsaids'}, axis=1, inplace=True)

    # mask for unspecified HRT answers
    HRT_mask = q142_cleaned.str.contains('hrt|estrogen|hormone replacement therapy|contracept')
    HRT_idx = q142_cleaned[HRT_mask].index.get_level_values(0)

    # mask for unspecified vitamin d3
    d3_mask = q142_cleaned.str.contains('(\s|^)d3|vitamin d(\s|$)', case=False)
    d3_idx = q142_cleaned[d3_mask].index.get_level_values(0)

    # mask for unspecified statin answers
    statin_mask = q142_cleaned.str.contains('(\s|^)statin(s|\s|$)', case=False)
    statin_idx = q142_cleaned[statin_mask].index.get_level_values(0)

    # mask for unspecified steroid answers
    steroid_mask = q142_cleaned.str.contains('(\s|^)corticosteroid(s|\s|$)', case=False)
    steroid_idx = q142_cleaned[steroid_mask].index.get_level_values(0)

    # set unspecified medications to NA
    patient_dose_feature_df.loc[HRT_idx, 'sex hormone therapy'] = 'NA'
    patient_dose_feature_df.loc[d3_idx, 'vitamin d and analogues'] = 'NA'
    patient_dose_feature_df.loc[statin_idx, 'statins'] = 'NA'
    patient_dose_feature_df.loc[steroid_idx, 'corticosteroids'] = 'NA'

    # save to csv file
    patient_dose_feature_df.to_csv('data/Covidence_12Aug20_Drug_Dosages.csv')