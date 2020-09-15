import pandas as pd
import warnings
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

# get patient numbers from the q142 cleaned object
patients = q142_cleaned.index.get_level_values(0).unique()

# drug classes and specific drugs to investigate
drug_classes = ['statins', 'ace inhibitors', 'proton pump inhibitors', 'corticosteroids', 'angiotensin ii receptor antagonists',
                'vitamin k antagonists', 'beta blocking agents', 'thiazides', 'h2-receptor antagonists', 'calcium-channel blockers', 'beta2-agonists',
                'antimuscarinics, other', 'non-steroidal anti-inflammatory drugs', 'sodium glucose co-transporter 2 inhibitors',
                'antiplatelet drugs', 'oestrogens|androgens', 'vitamin d and analogues', '^calcium$', 'bisphosphonates']

specific_drugs = ['paracetamol', 'metformin', 'aspirin']

# dictionary for patient drug classes
patient_drug_class_dict = {drug_class: [] for drug_class in drug_classes}

# dictionary for patient drugs
patient_drug_dict = {drug: [] for drug in specific_drugs}

# function for updating dictionary with a list of patient indices belonging to a drug class
def label_patient_classes(drug_class):

    # get bnf entries that fall into the desired drug class
    drug_class_mask = bnf_classes['primary'].str.contains(drug_class) | bnf_classes['secondary'].str.contains(drug_class)
    class_drugs = bnf_classes[drug_class_mask]

    # isolate entries that in the drug class that correspond to single drugbank IDs, to avoid ambiguity with mixture products
    single_id_mask = class_drugs['db_id'].apply(lambda ids: len(ids) == 1)
    # get the db ids corresponding to the drug class
    class_db_ids = set().union(*class_drugs['db_id'][single_id_mask])

    # loop through patients and check if any of the survey answers are in the set of class drugbank ids
    for patient in patients:
        patient_drugs = q142_cleaned[patient]
        is_in_class = patient_drugs.apply(lambda answer: any([db_id in class_db_ids for db_id in drug_dictionary.get(answer)]) if drug_dictionary.get(answer) else False)
        if is_in_class.any():
            patient_drug_class_dict[drug_class].append(patient)

# function for getting a list of patient indices taking a specific drug
def label_patient_drugs(drug):

    # get the drugbank id corresponding to the drug of interest
    drug_db_ids = drug_dictionary[drug]
    synonyms = set([name for name, db_ids in drug_dictionary.items() if drug_db_ids.issubset(db_ids)])

    # loop through patients and determine whether each person takes the drug
    for patient in patients:
        patient_drugs = q142_cleaned[patient]
        takes_drug = patient_drugs.apply(lambda drug: drug in synonyms).any()
        if takes_drug:
            patient_drug_dict[drug].append(patient)

# label patient drug classes
for drug_class in drug_classes:
    label_patient_classes(drug_class)

# label patient drugs
for drug in specific_drugs:
    label_patient_drugs(drug)

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
d3_mask = q142_cleaned.str.contains('(\s|^)d3|vitamin d|evacal (d3)?(\s|$)', case = False)
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

patient_feature_df.to_csv('data/Covidence_12Aug20_Drug_Classes.csv', index = False)