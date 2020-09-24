import pandas as pd
import re
import pickle
from abydos.phonetic import Metaphone

# import the drug dictionary
drug_dictionary = pickle.load(open('data/drug_dictionary.p', 'rb'))

# import the survey data
survey_data = pd.read_csv('data/Covidence_12Aug20_DrgExtra.csv')

# filter for medication question
q142 = survey_data.loc[:, survey_data.columns.str.contains('q142')]
# mask for rows (participants) that provided no medication answers
no_q142_answer_mask = q142.apply(lambda row: row.isna().all(), axis=1)
# filter participants that provided no answer and flatten the array
q142_filtered = q142[~no_q142_answer_mask].stack()
q142_filtered = q142_filtered[(q142_filtered != '-99') & (q142_filtered != -99)]
# basic preprocessing
q142_filtered = q142_filtered.str.strip()
q142_filtered = q142_filtered.str.lower()

# filter for dosage question
q143 = survey_data.loc[:, survey_data.columns.str.contains('q143')]
no_q143_answer_mask = q143.apply(lambda row: row.isna().all(), axis=1)
q143_filtered = q143[~no_q143_answer_mask].stack()

# separate dosage values and units
q143_dosages = q143_filtered[q143_filtered.index.get_level_values(1).str.contains('q1431')]
q143_units = q143_filtered[q143_filtered.index.get_level_values(1).str.contains('q1432')]

# rename the question index in the q143_units question to match the medication/dosage questions
q143_units.rename({q_id: q_id + '_1' for q_id in q143_units.index.get_level_values(1).unique()}, level = 1, inplace = True)

# add dosage answers for patients that provided medication answers but are missing in the dosage answers
for patient, question in q142_filtered.index:
    # get question index values for the dosage/q143_units questions
    dosage_question = question.replace('q1421', 'q1431')
    units_question = question.replace('q1421', 'q1432')
    # add missing patient/question combinations to the series
    if (patient, dosage_question) not in q143_dosages.index:
        q143_dosages[(patient, dosage_question)] = -99
        q143_dosages = q143_dosages.sort_index()
    if (patient, units_question) not in q143_units.index:
        q143_units[(patient, units_question)] = -99
        q143_units = q143_units.sort_index()

# generator to avoid modifying the global object
dosage_generator = ((patient, dosage_question) for patient, dosage_question in q143_dosages.index)
# for patients/questions that are in the dosage answer series but not the medication answer series, just drop them
for patient, dosage_question in dosage_generator:
    question = dosage_question.replace('q1431', 'q1421')
    units_question = dosage_question.replace('q1431', 'q1432')
    if (patient, question) not in q142_filtered.index:
        q143_dosages = q143_dosages.drop((patient, dosage_question))
    if (patient, units_question) not in q143_units.index:
        q143_dosages = q143_dosages.drop((patient, dosage_question))

# generator to avoid modifying the global object
units_generator = ((patient, units_question) for patient, units_question in q143_units.index)
# same for patients/questions in the dosage unit question
for patient, units_question in units_generator:
    question = units_question.replace('q1432', 'q1421')
    dosage_question = units_question.replace('q1432', 'q1431')
    if (patient, question) not in q142_filtered.index:
        q143_units = q143_units.drop((patient, units_question))
    if (patient, dosage_question) not in q143_dosages.index:
        q143_units = q143_units.drop((patient, units_question))

## regex patterns for cleaning medication answers ##

# pattern for drug weights (e.g. milligrams)
weights = 'm?(milli)?(micro)?(mc)?(mic)?\s*g(ram)?'
# pattern for volumes
volumes = 'm?(milli)?(micro)?(mc)?(mic)?\s*l(iter)?'
# weight and volume patterns combined
units = '(({weight}|{volume}|%|unit|i\.*u\.*)s*)'.format(weight=weights, volume=volumes)
# q143_units pattern compiled into a regex pattern for q143_dosages
dosage_pattern = '([\d.]+/)*([\d.]+\s*|\s+){units}((/|{units})|\s+|$)|(\s+[\d.x]+\s*/\s*[\d.]+(\s+|$))'.format(units=units)
dosage_regex = re.compile(dosage_pattern)

# pattern for different drug formulations
formulation_pattern = '(^|\s+)(capsule|drop|cream|ointment|tab(let)*|lotion|pill|spray|shampoo|patch(e)*|inhaler|gel|injection|pump|pen|solution|aqueous|oil|app(lication)*|implant|foam)s*(\s+|$)'
formulation_regex = re.compile(formulation_pattern)

# pattern for different routes of administration
routes_of_admin_pattern = '(^|\s+)((oral|nasal|ocular|auricular|topical)(ly)?|mouth|eye|nose|ear|skin|scalp)(\s+|$)'
RoA_regex = re.compile(routes_of_admin_pattern)

# pattern for numbers
numbers_pattern = '(once|one|1)|(twice|two|2)|(three|3)|(four|4)|(five|5)|(six|6)|(seven|7)|(eight|8)|(nine|9)'
# pattern for frequency of taking drugs (e.g. 2x a day)
frequency_pattern = '({numbers})?\s*(times|x)?\s*(a|per|every|each)?\s*({numbers})?\s*(da(y|ily)|(week|month)(ly)?)'.format(numbers=numbers_pattern)
frequency_regex = re.compile(frequency_pattern)

# pattern for parenthetical qualifiers
qualifier_regex = re.compile('\(+.*\)+')

all_patterns = [frequency_regex, dosage_regex, formulation_regex, RoA_regex, qualifier_regex]

# loop through regex patterns and filter each one from the answers
for pattern in all_patterns:
    q142_filtered = q142_filtered.str.replace(pattern, '')

# remove any thing coming after a forward slash if more than two alphanumeric characters are detected
q142_cleaned = q142_filtered.str.replace('/\w{2,}.*$', '')

# post processing
q142_cleaned = q142_cleaned.str.strip()
q142_cleaned = q142_cleaned.str.lower()

## make objects for tracking answer mappings

# objects for tracking the frequency of drugbank ids in the drug dictionary
all_db_ids = set().union(*drug_dictionary.values())
drug_frequencies = {db_id: 0 for db_id in all_db_ids}

# lists for mapped answers in different categories
mapped_survey_answers = []
first_name_mapped_survey_answers = []
unmapped_survey_answers = []

# loop through answers and map them
for answer in q142_cleaned:

    # try to get the drugbank ids for the whole answer
    db_ids = drug_dictionary.get(answer)

    # regex pattern to isolate first word
    first_word = re.sub('[^\w]+.*$', '', answer)
    first_word_db_ids = drug_dictionary.get(first_word)

    # if the name is already in the drug dictionary add to the mapped list
    if db_ids:
        mapped_survey_answers.append(answer)
        mapped_db_ids = db_ids

    # if its first name is in the bnf add it to the first name mapped list and add answer to the drug dictionary
    elif first_word_db_ids:
        drug_dictionary[answer] = drug_dictionary[first_word]
        first_name_mapped_survey_answers.append(answer)
        mapped_db_ids = first_word_db_ids

    # otherwise add it to the unmapped list
    else:
        unmapped_survey_answers.append(answer)
        mapped_db_ids = set()
    # for each of the drugbank ids, update the frequency dictionary
    for db_id in mapped_db_ids:
        drug_frequencies[db_id] += 1

## use metaphone to map phonetic encodings to drugbank ids in the drug dictionaries ##

mp = Metaphone()
# dictionary for storing drugs mapped to their encodings
drug_encodings = {}
# dictionary for storing encodings mapped to drugbank ids
encoded_drug_dict = {}
# list for ambigious encodings (distinct phonetically-identical drugs) - these will be removed from the dictionary
ambiguous_encodings = []

# loop through the drug dictionary and encode every entry, saving the corresponding drugbank ids under the encoding
for drug in drug_dictionary:

    # save the encoding for each drug
    encoding = mp.encode(drug)
    drug_encodings[drug] = encoding

    # if the encoding is not in the encoding dictionary, add it
    if encoding not in encoded_drug_dict:
        encoded_drug_dict[encoding] = drug_dictionary[drug]

    # if the encoding is already in the dictionary and there exists different ids for the same encoding, save it
    elif drug_dictionary[drug] != encoded_drug_dict[encoding]:
        ambiguous_encodings.append(encoding)

# filter for encodings with only one match in the drugbank
encoded_drug_dict = {key: val for key, val in encoded_drug_dict.items() if key not in ambiguous_encodings}

# get survey answers whose encodings are valid
mapped_by_encoding = [answer for answer in unmapped_survey_answers if mp.encode(answer) in encoded_drug_dict]
# add answers to the drug dictionary under the encoding's drugbank ids
for answer in mapped_by_encoding:
    drug_dictionary[answer] = encoded_drug_dict[mp.encode(answer)]

# list for drugs still unmapped by phonetic encoding
unmapped_by_encoding = [answer for answer in unmapped_survey_answers if answer not in mapped_by_encoding]

print('{} survey answers mapped'.format(len(mapped_survey_answers) + len(first_name_mapped_survey_answers) + len(mapped_by_encoding)))
print('{} survey answers unmapped'.format(len(unmapped_by_encoding)))

# get counts for different mapping categories
mapping_counts = {'exact': len(mapped_survey_answers), 'first_name': len(first_name_mapped_survey_answers),
                  'phonetic_encoding': len(mapped_by_encoding), 'unmapped': len(unmapped_by_encoding)}