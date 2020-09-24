from abydos.distance import Levenshtein
import numpy as np
import pandas as pd
import pickle

if __name__ == '__main__':

    # import unmapped answers
    from Map_survey_answers import AnswerMapper

    # load in drug dictionary
    drug_dictionary = pickle.load(open('data/drug_dictionary.p', 'rb'))

    # create instance of answer mapper class with the right survey file path
    mapper = AnswerMapper(drug_dictionary, survey_filepath='data/Covidence_12Aug20_DrgExtra.csv')

    # call map answers
    mapper.map_answers()

    drug_dictionary = mapper.drug_dictionary

    # final round of mapping, getting all drugs in the drug dictionary that are levenshtein distance of 1 from each answer
    lv = Levenshtein(mode = 'osa')
    mapped_by_lv_distance = {}
    unmapped_by_lv_distance = []
    # loop through unmapped answers
    for i, answer in enumerate(set(mapper.unmapped_by_encoding)):
        # only check for answers greater than 3 letters
        if len(answer) > 3:
            # do not map answers that are longer than one word
            if len(answer.split(' ')) > 1:
                unmapped_by_lv_distance.append(answer)
                continue
            # narrow down potential matches as those of similar length (plus or minus 3 characters) to speed up computation
            longest_len = len(answer)+1
            shortest_len = len(answer)-1
            potential_matches = (alias for alias in drug_dictionary if len(alias) >= shortest_len and len(alias) <= longest_len)
            # get distance with all potential matches
            distances = {alias: lv.dist_abs(answer, alias) for alias in potential_matches}
            # get matches that are a distance of 1 away
            matches = [alias for alias, distance in distances.items() if distance == 1]
            # if there are multiple matches that are a distance of 1 away, take the one that appears at the highest frequency
            if matches:
                best_match = max(matches, key = lambda alias: np.mean([drug_frequencies[db_id] for db_id in drug_dictionary[alias]]))
                mapped_by_lv_distance[answer] = best_match
            # if there are no matches 1 character away, save to the list
            else:
                unmapped_by_lv_distance.append(answer)
        else:
            unmapped_by_lv_distance.append(answer)
        # print every 100 answers, for checking progress
        if i % 100 == 0:
            print('Answer number {} completed'.format(i))

    # dump data frame of unmapped answers for manual annotation
    mapped_answer_df = pd.DataFrame([(answer, mapping) for answer, mapping in mapped_by_lv_distance.items()], columns = ['answer', 'correction'])
    unmapped_answer_df = pd.DataFrame([(answer, 0) for answer in unmapped_by_lv_distance], columns = ['answer', 'correction'])
    full_mapping_df = pd.concat([mapped_answer_df, unmapped_answer_df])
    full_mapping_df.to_csv('data/answer_mappings.csv', index = False)