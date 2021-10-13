import pandas as pd
import numpy as np
import warnings
import pickle
import re
import argparse
from itertools import compress
warnings.simplefilter('ignore')

# import class for mapping survey answers
from utils.answer_mapping import AnswerMapper

# import the drug dictionary
drug_dictionary = pickle.load(open('../data/drug_dictionary.p', 'rb'))

# drug classes and specific drugs to investigate
drug_classes = ['statins', 'ace inhibitors', 'proton pump inhibitors', 'corticosteroids',
                'angiotensin ii receptor antagonists', 'selective serotonin re-uptake inhibitors',
                'vitamin k antagonists', 'beta blocking agents', 'thiazides', 'h2-receptor antagonists',
                'calcium-channel blockers', 'beta2-agonists', 'antimetabolites',
                'calcineurin inhibitors', 'tumor necrosis factor alpha \(tnf-a\) inhibitors',
                'tricyclic antidepressants', 'serotonin and noradrenaline re-uptake inhibitors',
                'monoamine-oxidase a and b inhibitors, irreversible',
                'interleukin inhibitors', 'disease-modifying anti-rheumatic drugs',
                'antimuscarinics, other', 'non-steroidal anti-inflammatory drugs',
                'sodium glucose co-transporter 2 inhibitors',
                'antiplatelet drugs', 'oestrogens|androgens', 'vitamin d and analogues', '^calcium$',
                'bisphosphonates']

specific_drugs = ['paracetamol', 'metformin', 'aspirin', 'digoxin']

# class for annotating patients with BNF drug classes
class PatientAnnotator:

    # function for reading in BNF data and annotating table with DB ids
    def read_in_bnf(self, drug_dictionary):

        # import bnf class dataframe
        bnf_classes = pd.read_csv('data/bnf_drug_classifications.csv')
        bnf_classes['drugs'] = bnf_classes['drugs'].str.split('; ')

        # map bnf columns to drugbank ids
        bnf_classes['db_id'] = bnf_classes['drugs'].apply(
            lambda drugs: [drug_dictionary.get(drug) for drug in drugs if drug_dictionary.get(drug)])
        # take the union of the set of ids
        bnf_classes['db_id'] = bnf_classes['db_id'].apply(lambda ids: set().union(*ids))

        # set class attribute
        self.bnf_classes = bnf_classes

        # list of drugbank ids mapped to bnf classes
        self.bnf_db_ids = set().union(*self.bnf_classes['db_id'].tolist())

    # method for counting how many BNF entries were mapped to the drug dictionary
    def count_BNF_mappings(self):

        # test whether each bnf drug is in the drug dictionary
        not_found_in_drugbank = self.bnf_classes['drugs'].apply(
            lambda drugs: [drug not in self.drug_dictionary for drug in drugs])

        # get mapped and unmapped bnf drugs
        unmapped_bnf_drugs = set()
        mapped_bnf_drugs = set()
        for drug, mask in zip(self.bnf_classes['drugs'], not_found_in_drugbank):
            unmapped = set(compress(drug, mask))
            mapped = set(compress(drug, [not unmapped for unmapped in mask]))
            unmapped_bnf_drugs = unmapped_bnf_drugs.union(unmapped)
            mapped_bnf_drugs = mapped_bnf_drugs.union(mapped)

        # map first word of each drug entry in the bnf
        first_name_unmapped = []
        first_name_mapped = []
        for drug in unmapped_bnf_drugs:
            first_part = drug.split(' ')[0]
            if first_part in self.drug_dictionary:
                self.drug_dictionary[drug] = self.drug_dictionary[first_part]
                first_name_mapped.append(drug)
            else:
                first_name_unmapped.append(drug)

        print('{} BNF drugs mapped'.format(len(mapped_bnf_drugs) + len(first_name_mapped)))
        print('{} BNF drugs unmapped'.format(len(first_name_unmapped)))

        self.bnf_mapped = mapped_bnf_drugs
        self.bnf_first_name_mapped = first_name_mapped
        self.bnf_unmapped = first_name_unmapped

    # requires a series of patient answers and a dictionary of drug aliases mapped to DB ids
    def __init__(self, meds, RoAs, drug_dict):
        self.meds = meds
        self.RoAs = RoAs
        self.drug_dictionary = drug_dict
        # get bnf data from the data directory
        self.read_in_bnf(drug_dict)

    # function for updating dictionary with a list of patient indices belonging to a drug class
    def get_patients_in_class(self, drug_class, roa = None):

        valid_bnf_classes = self.bnf_classes[self.bnf_classes['db_id'].apply(lambda ids: len(ids) > 0)]

        # get bnf entries that fall into the desired drug class
        drug_class_mask = (valid_bnf_classes['primary'].str.contains(drug_class) |
                           valid_bnf_classes['secondary'].str.contains(drug_class))
        class_drugs = valid_bnf_classes[drug_class_mask]

        # isolate entries that in the drug class that correspond to single drugbank IDs, to avoid ambiguity with mixture products
        single_id_mask = class_drugs['db_id'].apply(lambda ids: len(ids) == 1)
        # get the db ids corresponding to the drug class
        class_db_ids = set().union(*class_drugs['db_id'][single_id_mask])

        # function for checking if a patient's answers are in a drug class
        # if roa not specified take all members from the class
        if roa is None:
            is_in_class = lambda patient: patient.apply(
                lambda answer: any([db_id in class_db_ids for db_id in self.drug_dictionary.get(answer)])
                if self.drug_dictionary.get(answer) else False).any()
            # group by patient and apply the is_in_class function to each patient
            patient_class_mask = self.meds.groupby(level=0).aggregate(func=is_in_class)
        # otherwise match by roa as well
        else:
            is_in_class = lambda meds, roas: any([any([db_id in class_db_ids for db_id in self.drug_dictionary.get(patient_med)])
                                                  and patient_roa == roa if self.drug_dictionary.get(patient_med) else False
                                                  for patient_med, patient_roa in zip(meds, roas)])
            patient_class_mask = []
            for patient in self.meds.index.get_level_values(0):
                in_class = is_in_class(self.meds[patient], self.RoAs[patient])
                patient_class_mask.append(in_class)
                
            patient_class_mask = pd.Series(patient_class_mask, index = self.meds.index.get_level_values(0))

        # get the indices of patients in the class
        patients_in_class = patient_class_mask.index[patient_class_mask].tolist()

        return patients_in_class

    # function for getting a list of patient indices taking a specific drug
    def get_patients_on_drug(self, drug):

        # get the drugbank ids corresponding to the drug of interest
        drug_db_ids = self.drug_dictionary[drug]

        # function for assessing whether each patient takes the drug
        takes_drug = lambda patient: patient.apply(lambda drug: self.drug_dictionary.get(drug).issuperset(drug_db_ids)
                                                   if self.drug_dictionary.get(drug) else False).any()

        # group by patient and apply the takes_drug function to each patient
        patient_drug_mask = self.meds.groupby(level = 0).aggregate(func = takes_drug)

        # get the indices of patients taking the drug
        patients_on_drug = patient_drug_mask.index[patient_drug_mask].tolist()

        return patients_on_drug

# if run from the command line, output a CSV file with answer mappings
if __name__ == '__main__':

    # file path and question column name arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('filepath', type=str, help='Path to the medication survey answers file')
    parser.add_argument('-q', '--questions', default = ['q1421', 'q1431', 'q1432', 'q1442'], nargs = 4, type = str,
                        help = 'Column names for medication, dosage, unit, and route of administration questions')
    parser.add_argument('-id', '--patient_id', default='uid', type=str, help='Column name for unique patient identifiers')
    args = parser.parse_args()

    # get the filename prefix from the filepath, for output file
    filename = re.search('.+(?=_.*\.csv$)', args.filepath).group(0)

    # create instance of answer mapper class with the survey file path
    mapper = AnswerMapper(survey_filepath=args.filepath, drug_dict=drug_dictionary, meds_q = args.questions[0],
                          dosage_q = args.questions[1], units_q = args.questions[2], RoAs_q = args.questions[3])

    # call map answers
    mapper.map_answers()

    # update drug dictionary with manual corrections file
    mapper.update_drug_dictionary(manual_corrections_filepath='data/answer_mappings_complete.csv')

    # dictionary for patient drug classes
    patient_drug_class_dict = {}

    # dictionary for patient drugs
    patient_drug_dict = {}

    # initialise class instance
    annotator = PatientAnnotator(meds=mapper.meds_cleaned, RoAs=mapper.RoAs, drug_dict=mapper.drug_dictionary)

    # label patient drug classes
    for drug_class in drug_classes:
        patient_drug_class_dict[drug_class] = annotator.get_patients_in_class(drug_class)

    patient_drug_class_dict['inhaled_corticosteroids'] = annotator.get_patients_in_class('corticosteroids', roa = 2)
    patient_drug_class_dict['oral_corticosteroids'] = annotator.get_patients_in_class('corticosteroids', roa = 1)

    # label patient drugs
    for drug in specific_drugs:
        patient_drug_dict[drug] = annotator.get_patients_on_drug(drug)

    # put all the patient features into a single dictionary
    patient_feature_dict = {**patient_drug_class_dict, **patient_drug_dict}

    # make a data frame
    patient_feature_df = pd.DataFrame(0, index = mapper.survey_data.index, columns = patient_feature_dict)
    patient_feature_df.insert(0, args.patient_id, mapper.survey_data[args.patient_id])

    # set patient indices for each feature to 1 in the full patient feature data frame
    for feature in patient_feature_dict:
        patient_feature_df.loc[patient_feature_dict[feature], feature] = 1

    # rename columns
    patient_feature_df.rename({'^calcium$': 'calcium', 'oestrogens|androgens': 'sex hormone therapy',
                               'antimuscarinics, other': 'antimuscarinics', 'non-steroidal anti-inflammatory drugs': 'nsaids',
                               'tumor necrosis factor alpha \(tnf-a\) inhibitors': 'tnf-a inhibitors'},
                              axis = 1, inplace = True)

    # manual corrections for HRT answers not in the drug dictionary
    HRT_mask = mapper.meds_cleaned.str.contains('hrt|estrogen|hormone replacement therapy|contracept')
    HRT_idx = mapper.meds_cleaned[HRT_mask].index.get_level_values(0)
    patient_feature_df.loc[HRT_idx, 'sex hormone therapy'] = 1

    # manual corrections for vitamin d3 answers not in the drug dictionary
    d3_mask = mapper.meds_cleaned.str.contains('(\s|^)d3|vitamin d(\s|$)', case = False)
    d3_idx = mapper.meds_cleaned[d3_mask].index.get_level_values(0)
    patient_feature_df.loc[d3_idx, 'vitamin d and analogues'] = 1

    # manual corrections for statin answers not in the drug dictionary
    statin_mask = mapper.meds_cleaned.str.contains('(\s|^)statin(s|\s|$)', case = False)
    statin_idx = mapper.meds_cleaned[statin_mask].index.get_level_values(0)
    patient_feature_df.loc[statin_idx, 'statins'] = 1

    # manual corrections for steroid answers not in the drug dictionary
    steroid_mask = mapper.meds_cleaned.str.contains('(\s|^)corticosteroid(s|\s|$)', case = False)
    steroid_idx = mapper.meds_cleaned[steroid_mask].index.get_level_values(0)
    patient_feature_df.loc[steroid_idx, 'corticosteroids'] = 1

    patient_feature_df['Systemic immunosuppressants'] = np.where((patient_feature_df['tnf-a inhibitors'] == 1) |
                                                                 (patient_feature_df['disease-modifying anti-rheumatic drugs'] == 1) |
                                                                 (patient_feature_df['interleukin inhibitors'] == 1) |
                                                                 (patient_feature_df['calcineurin inhibitors'] == 1) |
                                                                 (patient_feature_df['antimetabolites'] == 1) |
                                                                 (patient_feature_df['oral_corticosteroids'] == 1), 1, 0)

    patient_feature_df['Non-SSRI antidepressants'] = np.where((patient_feature_df['tricyclic antidepressants'] == 1) |
                                                              (patient_feature_df['monoamine-oxidase a and b inhibitors, irreversible'] == 1) |
                                                              (patient_feature_df['serotonin and noradrenaline re-uptake inhibitors'] == 1),
                                                              1, 0)

    # save the drug class data as a csv file
    patient_feature_df.to_csv('{}_Drug_Classes.csv'.format(filename), index = False)