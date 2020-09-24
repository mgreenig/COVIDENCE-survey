import pandas as pd
import re
from abydos.phonetic import Metaphone

# class for mapping survey answers
class AnswerMapper:

    # clean the survey answers
    def import_and_clean(self, survey_filepath):

        # import the survey data
        self.survey_data = pd.read_csv(survey_filepath)

        # filter for medication question
        meds = self.survey_data.loc[:, self.survey_data.columns.str.contains('q142')]
        # mask for rows (participants) that provided no medication answers
        no_q142_answer_mask = meds.apply(lambda row: row.isna().all(), axis=1)

        # filter participants that provided no answer and flatten the array
        meds_filtered = meds[~no_q142_answer_mask].stack()
        # filter for NAs or 99s
        meds_filtered = meds_filtered[(meds_filtered != '-99') & (meds_filtered != -99) & (~meds_filtered.isna())]
        # basic preprocessing
        meds_filtered = meds_filtered.str.strip()
        meds_filtered = meds_filtered.str.lower()

        # filter for dosage question
        dosages_and_units = self.survey_data.loc[:, self.survey_data.columns.str.contains('q143')]
        no_dosage_unit_answer_mask = dosages_and_units.apply(lambda row: row.isna().all(), axis=1)
        dosages_and_units_filtered = dosages_and_units[~no_dosage_unit_answer_mask].stack()

        # separate dosage values and units
        dosages = dosages_and_units_filtered[dosages_and_units_filtered.index.get_level_values(1).str.contains('q1431')]
        units = dosages_and_units_filtered[dosages_and_units_filtered.index.get_level_values(1).str.contains('q1432')]

        # rename the question index in the q143_units question to match the medication/dosage questions
        units.rename({q_id: q_id + '_1' for q_id in units.index.get_level_values(1).unique()}, level = 1, inplace = True)

        # add dosage answers for patients that provided medication answers but are missing in the dosage answers
        for patient, question in meds_filtered.index:
            # get question index values for the dosage/q143_units questions
            dosage_question = question.replace('q1421', 'q1431')
            units_question = question.replace('q1421', 'q1432')
            # add missing patient/question combinations to the series
            if (patient, dosage_question) not in dosages.index:
                dosages[(patient, dosage_question)] = -99
                dosages = dosages.sort_index()
            if (patient, units_question) not in units.index:
                units[(patient, units_question)] = -99
                units = units.sort_index()

        # generator to avoid modifying the global object
        dosage_generator = ((patient, dosage_question) for patient, dosage_question in dosages.index)
        # for patients/questions that are in the dosage answer series but not the medication answer series, just drop them
        for patient, dosage_question in dosage_generator:
            question = dosage_question.replace('q1431', 'q1421')
            units_question = dosage_question.replace('q1431', 'q1432')
            if (patient, question) not in meds_filtered.index:
                dosages = dosages.drop((patient, dosage_question))
            if (patient, units_question) not in units.index:
                dosages = dosages.drop((patient, dosage_question))

        # generator to avoid modifying the global object
        units_generator = ((patient, units_question) for patient, units_question in units.index)
        # same for patients/questions in the dosage unit question
        for patient, units_question in units_generator:
            question = units_question.replace('q1432', 'q1421')
            dosage_question = units_question.replace('q1432', 'q1431')
            if (patient, question) not in meds_filtered.index:
                units = units.drop((patient, units_question))
            if (patient, dosage_question) not in dosages.index:
                units = units.drop((patient, units_question))

        ## regex patterns for cleaning medication answers ##

        # pattern for drug weights (e.g. milligrams)
        weights = 'm?(milli)?(micro)?(mc)?(mic)?\s*g(ram)?'
        # pattern for volumes
        volumes = 'm?(milli)?(micro)?(mc)?(mic)?\s*l(iter)?'
        # weight and volume patterns combined
        dose_units = '(({weight}|{volume}|%|unit|i\.*u\.*)s*)'.format(weight=weights, volume=volumes)
        # q143_units pattern compiled into a regex pattern for q143_dosages
        dosage_pattern = '([\d.]+/)*([\d.]+\s*|\s+){units}((/|{units})|\s+|$)|(\s+[\d.x]+\s*/\s*[\d.]+(\s+|$))'.format(units=dose_units)
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
            meds_filtered = meds_filtered.str.replace(pattern, '')

        # remove any thing coming after a forward slash if more than two alphanumeric characters are detected
        meds_cleaned = meds_filtered.str.replace('/\w{2,}.*$', '')

        # post processing
        meds_cleaned = meds_cleaned.str.strip()
        meds_cleaned = meds_cleaned.str.lower()

        self.meds_cleaned = meds_cleaned
        self.dosages = dosages
        self.units = units

    # initialise with a drug dictionary and a filepath to the survey data frame
    def __init__(self, survey_filepath, drug_dict):
        self.drug_dictionary = drug_dict
        self.all_db_ids = set().union(*self.drug_dictionary.values())
        self.drug_frequencies = {db_id: 0 for db_id in self.all_db_ids}
        self.import_and_clean(survey_filepath)

    def map_answers(self):

        # lists for mapped answers in different categories
        self.mapped_survey_answers = []
        self.first_name_mapped_survey_answers = []
        self.unmapped_survey_answers = []

        # loop through answers and map them
        for answer in self.meds_cleaned:

            # try to get the drugbank ids for the whole answer
            db_ids = self.drug_dictionary.get(answer)

            # regex pattern to isolate first word
            first_word = re.sub('[^\w]+.*$', '', answer)
            first_word_db_ids = self.drug_dictionary.get(first_word)

            # if the name is already in the drug dictionary add to the mapped list
            if db_ids:
                self.mapped_survey_answers.append(answer)
                mapped_db_ids = db_ids

            # if its first name is in the bnf add it to the first name mapped list and add answer to the drug dictionary
            elif first_word_db_ids:
                self.drug_dictionary[answer] = self.drug_dictionary[first_word]
                self.first_name_mapped_survey_answers.append(answer)
                mapped_db_ids = first_word_db_ids

            # otherwise add it to the unmapped list
            else:
                self.unmapped_survey_answers.append(answer)
                mapped_db_ids = set()
            # for each of the drugbank ids, update the frequency dictionary
            for db_id in mapped_db_ids:
                self.drug_frequencies[db_id] += 1

        ## use metaphone to map phonetic encodings to drugbank ids in the drug dictionaries ##

        mp = Metaphone()
        # dictionary for storing encodings mapped to drugbank ids
        encoded_drug_dict = {}
        # list for ambigious encodings (distinct phonetically-identical drugs) - these will be removed from the dictionary
        ambiguous_encodings = []

        # loop through the drug dictionary and encode every entry, saving the corresponding drugbank ids under the encoding
        for drug in self.drug_dictionary:

            # save the encoding for each drug
            encoding = mp.encode(drug)

            # if the encoding is not in the encoding dictionary, add it
            if encoding not in encoded_drug_dict:
                encoded_drug_dict[encoding] = self.drug_dictionary[drug]

            # if the encoding is already in the dictionary and there exists different ids for the same encoding, save it
            elif self.drug_dictionary[drug] != encoded_drug_dict[encoding]:
                ambiguous_encodings.append(encoding)

        # filter for encodings with only one match in the drugbank
        encoded_drug_dict = {key: val for key, val in encoded_drug_dict.items() if key not in ambiguous_encodings}

        # get survey answers whose encodings are valid
        self.mapped_by_encoding = [answer for answer in self.unmapped_survey_answers if mp.encode(answer) in encoded_drug_dict]
        # add answers to the drug dictionary under the encoding's drugbank ids
        for answer in self.mapped_by_encoding:
            self.drug_dictionary[answer] = encoded_drug_dict[mp.encode(answer)]

        # list for drugs still unmapped by phonetic encoding
        self.unmapped_by_encoding = [answer for answer in self.unmapped_survey_answers if answer not in self.mapped_by_encoding]

    def update_drug_dictionary(self, manual_corrections_filepath):

        # manual id setting for a few medications
        self.drug_dictionary['vitamin b12'] = {'DB00115'}
        self.drug_dictionary['vitamin e'] = {'DB00163'}
        self.drug_dictionary['vitamin d'] = {'DB00136', 'DB00153', 'DB00169', 'DB00910', 'DB01070', 'DB01436', 'DB02300',
                                        'DB13689'}
        self.drug_dictionary['candesartan cilexetil'] = {'DB13919'}
        self.drug_dictionary['candesartan'] = {'DB13919'}

        # load in unmapped answers and map to the drug dictionary
        corrections = pd.read_csv(manual_corrections_filepath)
        corrections = corrections.astype(str)
        corrections = corrections.apply(lambda col: col.str.strip(), axis=0)

        # add unmapped answers to drug dictionary
        for answer, correction in zip(corrections['answer'], corrections['correction']):
            correction_split = correction.split('; ')
            if any([corr in self.drug_dictionary for corr in correction_split]):
                self.drug_dictionary[answer] = set().union(
                    *[self.drug_dictionary.get(corr) for corr in correction_split if self.drug_dictionary.get(corr)])
            elif str(correction) != '0':
                self.meds_cleaned[self.meds_cleaned == answer] = correction