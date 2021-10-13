import pandas as pd
import numpy as np
import re
import argparse
import pickle
from scipy.stats import zscore, norm

# import the drug dictionary
drug_dictionary = pickle.load(open('../data/drug_dictionary.p', 'rb'))

# import the relevant objects
from Map_survey_answers import AnswerMapper
from Annotate_patients import PatientAnnotator, drug_classes, specific_drugs

# class for scaling dosage values
class DosageScaler(PatientAnnotator):

    # class takes in the full survey answers, medication answers, dosage answers, dosage units answers, and a drug dictionary
    def __init__(self, survey_data, meds, dosages, units, drug_dict):
        '''
        The meds, dosages, and units series should be on the same multi-index of the form (patient_number, question_number)

        We use the first five letters of the question number index to distinguish between med, dosage, and unit questions, i.e.
        - q1421_x -> meds
        - q1431_x -> dosage values
        - q1432_x -> dosage units
        '''
        # initialise parent class to call read_bnf()
        super().__init__(meds, drug_dict)
        self.survey_data = survey_data
        self.med_db_ids = meds.apply(lambda answer: drug_dict.get(answer) if drug_dict.get(answer) else set())
        self.dosages = dosages
        self.units = units

    # function for aligning question indices in a survey answer mask to another survey answer series
    @staticmethod
    def align_mask(mask, series):
        # index pattern is everything before the first underscore in the first index entry
        mask_idx_pattern = re.sub('_.*$', '', mask.index[0][1])
        series_idx_pattern = re.sub('_.*$', '', series.index[0][1])
        # get the question patterns from of the first five letters of the mask and answer series indices
        q1_pattern = mask.index[0][1][:min(5, len(mask_idx_pattern))]
        q2_pattern = series.index[0][1][:min(5, len(series_idx_pattern))]
        # rename the question index in the mask
        mask_renamed = mask.rename({q_id: q_id.replace(q1_pattern, q2_pattern) for q_id in mask.index.get_level_values(1).unique()}, level = 1)
        # align it to the series
        _, mask_aligned = series.align(mask_renamed, fill_value = False)
        return mask_aligned

    # function for getting normalised drug dosages for a drugbank ID
    def get_normalised_dosages(self, id):

        # if a single id is provided, convert into a set
        ids = set([id]) if isinstance(id, str) else id

        # make copy of the dosage question so the global object is not modified
        dosages = self.dosages.copy()
        # filter for answers containing the same DB id(s)
        drug_mask = self.med_db_ids.apply(lambda drug_ids: True if ids.issubset(drug_ids) else False)
        # change index name
        drug_mask = DosageScaler.align_mask(drug_mask, self.units)

        # if any drugs are hit, get the dosages
        if drug_mask.any():

            # filter for answers containing the same DB id(s)
            exact_drug_mask = self.med_db_ids.apply(lambda drug_ids: True if ids == drug_ids else False)
            # change index name
            exact_drug_mask = DosageScaler.align_mask(exact_drug_mask, self.units)

            # get mixture compounds
            mixture_mask = DosageScaler.align_mask(drug_mask != exact_drug_mask, dosages)

            # if there are any exact matches, scale the drug dosages
            if exact_drug_mask.any():

                # get all dosage units specified for the drug
                dosage_units = self.units[exact_drug_mask].copy()
                mg_mask = DosageScaler.align_mask(dosage_units == 1, dosages)
                micg_mask = DosageScaler.align_mask(dosage_units == 2, dosages)
                # get invalid units
                invalid_unit_mask = DosageScaler.align_mask((dosage_units != 1) & (dosage_units != 2), dosages)

                # get iqr and quantiles of the mg drug dosages
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
            else:
                invalid_unit_mask = pd.Series(False, index = dosages.index)
                valid_dosages_norm = pd.Series()

            # NA values are either mixtures or invalid dosages
            NA_mask = invalid_unit_mask | mixture_mask

            # get invalid dosages
            invalid_dosages = dosages[NA_mask].apply(lambda dosage: -1)

            # combine the two dosage lists
            all_dosage_values = invalid_dosages.append(valid_dosages_norm)
            all_dosage_values = all_dosage_values.sort_index()
            all_dosage_values[(all_dosage_values == -99) | (all_dosage_values == '-99')] = -1

            # align to the larger series
            _, drug_dosages_aligned = dosages.align(all_dosage_values, fill_value=0)

        # otherwise save the dosages as a series of 0
        else:
            drug_dosages_aligned = pd.Series(0, index = dosages.index)

        return drug_dosages_aligned

    @staticmethod
    # function for checking if a series of dosage answers consists of only NA and 0 values
    def is_na(row):
        return False if all(row == 0) else False if all(row != -1) else True

    @staticmethod
    # function for combining dosage data across multiple drugs into a single value (either a sum or NA)
    def combine_func(row):
        return -1 if DosageScaler.is_na(row) else row[row != -1].sum()

    # function for getting the dosage values for each class
    def get_class_doses(self, drug_class):

        # filter for BNF listings that mapped to the drug dictionary
        valid_class_mask = self.bnf_classes['db_id'].apply(lambda ids: len(ids) > 0)
        valid_bnf_classes = self.bnf_classes[valid_class_mask]

        # get bnf entries that fall into the desired drug class
        drug_class_mask = valid_bnf_classes['primary'].str.contains(drug_class) | valid_bnf_classes['secondary'].str.contains(drug_class)
        class_drugs = valid_bnf_classes[drug_class_mask]

        # isolate entries that in the drug class that correspond to single drugbank IDs, to avoid ambiguity with mixture products
        single_id_mask = class_drugs['db_id'].apply(lambda ids: len(ids) == 1)
        # get the db ids corresponding to the drug class
        class_db_ids = set().union(*class_drugs['db_id'][single_id_mask])

        # loop through ids in the class and add dosages to the dictionary of dosage data
        drug_dosages = {}
        for id in class_db_ids:
            dosages = self.get_normalised_dosages(id)
            drug_dosages[id] = dosages

        # make data frame from the dosage data for all drugs in the class
        dosage_df = pd.DataFrame(drug_dosages)

        # group by patient, summing to get the total dose per patient within each class
        question_dosages = dosage_df.apply(DosageScaler.combine_func, axis = 1)
        patient_dosages = question_dosages.groupby(level = 0).apply(DosageScaler.combine_func)

        # align to the total set of survey answers
        _, patient_dosages_aligned = self.survey_data.align(patient_dosages, axis = 0, fill_value = 0)

        return patient_dosages_aligned

    # function for getting drug dosages
    def get_drug_doses(self, drug):

        # get the drugbank id corresponding to the drug of interest
        drug_db_ids = drug_dictionary[drug]
        # get the normalised dosages
        dosages = self.get_normalised_dosages(drug_db_ids)
        # group by patient, summing to get the total dose per patient within each class
        patient_dosages = dosages.groupby(level = 0).apply(self.combine_func)
        # align to the total set of survey answers
        _, patient_dosages_aligned = self.survey_data.align(patient_dosages, axis = 0, fill_value = 0)

        return patient_dosages_aligned


if __name__ == '__main__':

    # file path and column name arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('filepath', type=str, help='Path to the medication survey answers file')
    parser.add_argument('-q', '--questions', default = ['q1421', 'q1431', 'q1432', 'q1442'], nargs = 4, type = str,
                        help = 'Column names for medication, dosage, and unit questions')
    parser.add_argument('-id', '--patient_id', default='uid', type=str, help='Column name for unique patient identifiers')
    args = parser.parse_args()

    # get the filename prefix from the filepath, for output file name
    filename = re.search('.+(?=_.*\.csv$)', args.filepath).group(0)

    # create instance of answer mapper class with the survey file path
    mapper = AnswerMapper(survey_filepath=args.filepath, drug_dict=drug_dictionary, meds_q = args.questions[0],
                          dosage_q = args.questions[1], units_q = args.questions[2], RoAs_q = args.questions[3])

    # generate answer mappings
    mapper.map_answers()

    # update drug dictionary with manual corrections file
    mapper.update_drug_dictionary(manual_corrections_filepath='../data/answer_mappings_complete.csv')

    # dictionary for patient drug classes
    drug_class_doses = {}

    # dictionary for patient drugs
    specific_drug_doses = {}

    # make a class instance with the mapped answer data
    scaler = DosageScaler(survey_data = mapper.survey_data, meds = mapper.meds_cleaned,
                          dosages = mapper.dosages, units = mapper.units, drug_dict = mapper.drug_dictionary)

    # label patient drug classes
    for drug_class in drug_classes:
        drug_class_doses[drug_class] = scaler.get_class_doses(drug_class)

    # label patient drugs
    for drug in specific_drugs:
        specific_drug_doses[drug] = scaler.get_drug_doses(drug)

    # put all the patient features into a single dictionary
    patient_dose_feature_dict = {**drug_class_doses, **specific_drug_doses}

    # make a data frame
    patient_dose_feature_df = pd.DataFrame(patient_dose_feature_dict)
    patient_dose_feature_df.insert(0, args.patient_id, mapper.survey_data[args.patient_id])

    patient_dose_feature_df.rename({'^calcium$': 'calcium', 'oestrogens|androgens': 'sex hormone therapy',
                                    'antimuscarinics, other': 'antimuscarinics',
                                    'non-steroidal anti-inflammatory drugs': 'nsaids'}, axis=1, inplace=True)

    # mask for unspecified HRT answers
    HRT_mask = mapper.meds_cleaned.str.contains('hrt|estrogen|hormone replacement therapy|contracept')
    HRT_idx = mapper.meds_cleaned[HRT_mask].index.get_level_values(0)

    # mask for unspecified vitamin d3
    d3_mask = mapper.meds_cleaned.str.contains('(\s|^)d3|vitamin d(\s|$)', case=False)
    d3_idx = mapper.meds_cleaned[d3_mask].index.get_level_values(0)

    # mask for unspecified statin answers
    statin_mask = mapper.meds_cleaned.str.contains('(\s|^)statin(s|\s|$)', case=False)
    statin_idx = mapper.meds_cleaned[statin_mask].index.get_level_values(0)

    # mask for unspecified steroid answers
    steroid_mask = mapper.meds_cleaned.str.contains('(\s|^)corticosteroid(s|\s|$)', case=False)
    steroid_idx = mapper.meds_cleaned[steroid_mask].index.get_level_values(0)

    # set unspecified medications to NA
    patient_dose_feature_df.loc[HRT_idx, 'sex hormone therapy'] = -1
    patient_dose_feature_df.loc[d3_idx, 'vitamin d and analogues'] = -1
    patient_dose_feature_df.loc[statin_idx, 'statins'] = -1
    patient_dose_feature_df.loc[steroid_idx, 'corticosteroids'] = -1

    # save to csv file
    patient_dose_feature_df.to_csv('../{}_Drug_Dosages.csv'.format(filename), index = False)