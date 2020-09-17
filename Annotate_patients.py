import pandas as pd
import warnings
import numpy as np
warnings.simplefilter('ignore')

from Map_survey_answers import q142_cleaned, survey_data, drug_dictionary

# import bnf class dataframe
bnf_classes = pd.read_csv('data/bnf_drug_classifications.csv')
# filter for drugs with drugbank ids
bnf_classes = bnf_classes[~bnf_classes['db_id'].isna()]
# split the drugbank ids to produce sets
bnf_classes['db_id'] = bnf_classes['db_id'].apply(lambda ids: set(ids.split('; ')))
# list of drugbank ids mapped to bnf classes
bnf_db_ids = set().union(*bnf_classes['db_id'].tolist())

# manual id setting for vitamins
drug_dictionary['vitamin b12'] = {'DB00115'}
drug_dictionary['vitamin e'] = {'DB00163'}
drug_dictionary['vitamin d'] = {'DB00136', 'DB00153', 'DB00169', 'DB00910', 'DB01070', 'DB01436', 'DB02300', 'DB13689'}

# load in unmapped answers and map to the drug dictionary
corrections = pd.read_csv('data/answer_mappings_complete.csv')
corrections = corrections.astype(str)
corrections = corrections.apply(lambda col: col.str.strip(), axis = 0)

# add unmapped answers to drug dictionary
for answer, correction in zip(corrections['answer'], corrections['correction']):
    correction_split = correction.split('; ')
    if any([corr in drug_dictionary for corr in correction_split]):
        drug_dictionary[answer] = set().union(*[drug_dictionary.get(corr) for corr in correction_split if drug_dictionary.get(corr)])
    elif str(correction) != '0':
        q142_cleaned[q142_cleaned == answer] = correction

# drug classes and specific drugs to investigate
drug_classes = ['statins', 'ace inhibitors', 'proton pump inhibitors', 'corticosteroids', 'angiotensin ii receptor antagonists',
                'vitamin k antagonists', 'beta blocking agents', 'thiazides', 'h2-receptor antagonists', 'calcium-channel blockers', 'beta2-agonists',
                'antimuscarinics, other', 'non-steroidal anti-inflammatory drugs', 'sodium glucose co-transporter 2 inhibitors',
                'antiplatelet drugs', 'oestrogens|androgens', 'vitamin d and analogues', '^calcium$', 'bisphosphonates']

specific_drugs = ['paracetamol', 'metformin', 'aspirin']

# dictionary for patient drug classes
patient_drug_class_dict = {}

# dictionary for patient drugs
patient_drug_dict = {}

# function for updating dictionary with a list of patient indices belonging to a drug class
def get_patients_in_class(drug_class, patient_answers):

    # get bnf entries that fall into the desired drug class
    drug_class_mask = bnf_classes['primary'].str.contains(drug_class) | bnf_classes['secondary'].str.contains(drug_class)
    class_drugs = bnf_classes[drug_class_mask]

    # isolate entries that in the drug class that correspond to single drugbank IDs, to avoid ambiguity with mixture products
    single_id_mask = class_drugs['db_id'].apply(lambda ids: len(ids) == 1)
    # get the db ids corresponding to the drug class
    class_db_ids = set().union(*class_drugs['db_id'][single_id_mask])

    # function for checking if a patient's answers are in a drug class
    is_in_class = lambda patient: patient.apply(
        lambda answer: any([db_id in class_db_ids for db_id in drug_dictionary.get(answer)]) if drug_dictionary.get(
            answer) else False).any()

    # group by patient and apply the is_in_class function to each patient
    patient_class_mask = patient_answers.groupby(level = 0).aggregate(func = is_in_class)

    # get the indices of patients in the class
    patients_in_class = patient_class_mask.index[patient_class_mask].tolist()

    return patients_in_class

# function for getting a list of patient indices taking a specific drug
def get_patients_on_drug(drug, patient_answers):

    # get the drugbank id corresponding to the drug of interest
    drug_db_ids = drug_dictionary[drug]

    # function for assessing whether each patient takes the drug
    takes_drug = lambda patient: patient.apply(lambda drug: drug_dictionary.get(drug).issuperset(drug_db_ids)
                                               if drug_dictionary.get(drug) else False).any()

    # group by patient and apply the takes_drug function to each patient
    patient_drug_mask = patient_answers.groupby(level = 0).aggregate(func = takes_drug)

    # get the indices of patients taking the drug
    patients_on_drug = patient_drug_mask.index[patient_drug_mask].tolist()

    return patients_on_drug

if __name__ == '__main__':

    # label patient drug classes
    for drug_class in drug_classes:
        patient_drug_class_dict[drug_class] = get_patients_in_class(drug_class, q142_cleaned)

    # label patient drugs
    for drug in specific_drugs:
        patient_drug_dict[drug] = get_patients_on_drug(drug, q142_cleaned)

    # put all the patient features into a single dictionary
    patient_feature_dict = {**patient_drug_class_dict, **patient_drug_dict}

    # make a data frame
    patient_feature_df = pd.DataFrame(0, index = survey_data.index, columns = patient_feature_dict)
    patient_feature_df.insert(0, 'uid', survey_data['uid'])

    for feature in patient_feature_dict:
        patient_feature_df.loc[patient_feature_dict[feature], feature] = 1

    patient_feature_df.rename({'^calcium$': 'calcium', 'oestrogens|androgens': 'sex hormone therapy',
                               'antimuscarinics, other': 'antimuscarinics', 'non-steroidal anti-inflammatory drugs': 'nsaids'},
                              axis = 1, inplace = True)

    # manual corrections for HRT answers not in the drug dictionary
    HRT_mask = q142_cleaned.str.contains('hrt|estrogen|hormone replacement therapy|contracept')
    HRT_idx = q142_cleaned[HRT_mask].index.get_level_values(0)
    patient_feature_df.loc[HRT_idx, 'sex hormone therapy'] = 1

    # manual corrections for vitamin d3 answers not in the drug dictionary
    d3_mask = q142_cleaned.str.contains('(\s|^)d3|vitamin d(\s|$)', case = False)
    d3_idx = q142_cleaned[d3_mask].index.get_level_values(0)
    patient_feature_df.loc[d3_idx, 'vitamin d and analogues'] = 1

    # manual corrections for statin answers not in the drug dictionary
    statin_mask = q142_cleaned.str.contains('(\s|^)statin(s|\s|$)', case = False)
    statin_idx = q142_cleaned[statin_mask].index.get_level_values(0)
    patient_feature_df.loc[statin_idx, 'statins'] = 1

    # manual corrections for steroid answers not in the drug dictionary
    steroid_mask = q142_cleaned.str.contains('(\s|^)corticosteroid(s|\s|$)', case = False)
    steroid_idx = q142_cleaned[steroid_mask].index.get_level_values(0)
    patient_feature_df.loc[steroid_idx, 'corticosteroids'] = 1

    # save the drug class data as a csv file
    patient_feature_df.to_csv('data/Covidence_12Aug20_Drug_Classes.csv', index = False)