import pandas as pd
import pickle

# import data and convert all to lowercase string
corrections = pd.read_csv('data/survey_answer_manual_corrections.csv', dtype = 'object')
n_corrections = len(corrections)
corrections = corrections.astype(str)
corrections = corrections.apply(lambda col: col.str.lower(), axis = 1)
corrections[corrections['Best guess (1)'] == 'nan'] = '5'
corrections[corrections['Best guess (2)'] == 'nan'] = '5'

# take first word of all custom answers
corrections['Candidate 4 (1)'] = corrections['Candidate 4 (1)'].str.split(' ').str[0]
corrections['Candidate 4 (2)'] = corrections['Candidate 4 (2)'].str.split(' ').str[0]

# mask for rows in which the mapping was uncertain
no_answer_mask = (corrections['Best guess (1)'] == '5') & (corrections['Best guess (2)'] == '5')
corrections = corrections[~no_answer_mask]

# mask for rows with custom answers
custom_answer_mask = (corrections['Best guess (1)'] == '4') & (corrections['Best guess (2)'] == '4')

# mask to filter custom answers that do not match
non_matching_custom_answer_mask = corrections['Candidate 4 (1)'] != corrections['Candidate 4 (2)']

# combine the two masks to filter rows in which custom answers were provided but do not match
ambiguous_custom_answer_mask = custom_answer_mask & non_matching_custom_answer_mask
corrections = corrections[~ambiguous_custom_answer_mask]

# mask to filter for rows in which best guesses do not match
matching_answer_mask = corrections['Best guess (1)'] == corrections['Best guess (2)']

# mask to filter for rows in which one guess is NA but the other one is provided
one_guess_NA_mask = (corrections['Best guess (1)'] == '5') ^ (corrections['Best guess (2)'] == '5')

# combine the two masks to filter for rows in which either both guesses match or one of the two guesses is NA
valid_guess_mask = matching_answer_mask | one_guess_NA_mask
corrections = corrections[valid_guess_mask]

# import the drug dictionary
drug_dictionary = pickle.load(open('data/drug_dictionary.p', 'rb'))

# function for getting the answer/correction pair for each row (to be added to the drug dictionary)
def map_corrections(row):
    # if either correction is 5, take the correction that is not equal to 5
    if (row['Best guess (1)'] == '5') ^ (row['Best guess (2)'] == '5'):
        if row['Best guess (1)'] == '5':
            col = 'Candidate ' + row['Best guess (2)']
        else:
            col = 'Candidate ' + row['Best guess (1)']
    # other take the correction specified by 'Best guess'
    else:
        # first assert the best guesses are the same
        assert row['Best guess (1)'] == row['Best guess (2)']
        col = 'Candidate ' + row['Best guess (1)']
    # if the column is Candidate 4 the name must be changed to match the table
    if col == 'Candidate 4':
        # check which of the two guesses are custom
        custom_guess_mask = row[['Best guess (1)', 'Best guess (2)']] == '4'
        # if both are custom take the first one
        if custom_guess_mask.all():
            col = 'Candidate 4 (1)'
        # else if only one is custom, set that custom answer as the column to be taken
        else:
            col = [colname for colname, is_custom_guess in zip(['Candidate 4 (1)', 'Candidate 4 (2)'], custom_guess_mask) if is_custom_guess]
            col = col[0]

    return row['Answer'], row[col]

# get the answer/correction pairs
mapped_corrections = corrections.apply(map_corrections, axis = 1)
corrections_df = pd.DataFrame(mapped_corrections.tolist(), columns = ['answer', 'correction'])
corrections_df.drop_duplicates(inplace = True)
corrections_df.to_csv('data/survey_answer_corrections.csv')

# add answers to the drug dictionary
n_corrections_mapped = 0
for answer, correction in mapped_corrections:
    if drug_dictionary.get(correction):
        drug_dictionary[answer] = drug_dictionary[correction]
        n_corrections_mapped += 1

print('{}/{} misspellings mapped using the manual corrections'.format(n_corrections_mapped, n_corrections))

pickle.dump(drug_dictionary, open('data/drug_dictionary.p', 'wb'))