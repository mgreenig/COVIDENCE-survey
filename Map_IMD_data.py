import pandas as pd
from bs4 import BeautifulSoup as bs
from urllib.request import Request, urlopen
import numpy as np
from multiprocessing import Pool, cpu_count

# IMD deciles for Northern irish postcodes
NI_imd_deciles = np.array([np.percentile(np.linspace(0, 890, 890), i) for i in range(0, 100, 10)])

# function for getting the IMD decile from a rank
def get_decile(rank):
    if rank:
        decile = np.count_nonzero(NI_imd_deciles < int(rank))
        return int(decile+1)
    else:
        return None

# mapping northern ireland postcodes to IMDs - have to pull from the web
def get_html(url):
    req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    webpage = urlopen(req).read()
    html = bs(webpage, 'html.parser')
    return html

# url for the northern irish IMD postcode lookup tool
base_url = 'https://deprivation.nisra.gov.uk/MDM/Details?Id='

# function for getting the postcode rank of a northern-irish postcode
def get_postcode_rank(postcode):
    postcode_for_url = postcode.replace(' ', '+')
    url = base_url + postcode_for_url
    html = get_html(url)
    imd_rank = html.find('span', {'style': 'font-size:48px; font-weight:bold;'})
    return imd_rank

if __name__ == '__main__':

    # import postcode data set
    postcode_data = pd.read_csv('data/postcode_data.csv', usecols = ['Postcode', 'In Use?', 'Country'])
    postcodes = postcode_data['Postcode']

    # import postcodes from COVIDENCE survey
    survey_postcodes = pd.read_csv('data/Covidence_12Aug20_Postcodes.csv')
    # remove trailing whitespaces
    survey_postcodes['pcode'] = survey_postcodes['pcode'].str.strip()
    # remove punctuation from the postcodes
    survey_postcodes['pcode'] = survey_postcodes['pcode'].str.replace('[^\w\s]', '')

    # dictionary for storing postcodes without spaces
    postcodes_space_removed = {}
    for postcode in postcodes:
        # remove spaces from the postcode
        postcode_joined = postcode.replace(' ', '')
        # if the postcode without spaces is not in the dictionary, add it
        if postcode_joined not in postcodes_space_removed:
            postcodes_space_removed[postcode_joined] = [postcode]
        # otherwise append it
        else:
            postcodes_space_removed[postcode_joined].append(postcode)

    # isolate survey postcodes that do not map to the postcode data set
    unmapped_postcodes = survey_postcodes.loc[~survey_postcodes['pcode'].isin(postcodes), 'pcode']
    # remove spaces
    unmapped_postcodes_space_removed = unmapped_postcodes.str.replace(' ', '')

    # go through unmapped postcodes and map to the postcode dictionary
    postcode_mappings = {}
    for postcode in unmapped_postcodes_space_removed:
        if postcodes_space_removed.get(postcode):
            # only map if there is a single possibility for the space-removed postcode
            if len(postcodes_space_removed[postcode]) == 1:
                postcode_mappings[postcode] = postcodes_space_removed[postcode][0]

    # indices and masks for postcodes mappable with the dictionary
    mappable_pcode_mask = unmapped_postcodes_space_removed.apply(lambda postcode: postcode in postcode_mappings)
    mappable_pcode_idx = mappable_pcode_mask.index[mappable_pcode_mask].values
    unmappable_pcode_idx = mappable_pcode_mask.index[~mappable_pcode_mask].values

    # map the postcodes using the mapping dictionary
    mapped_postcodes = unmapped_postcodes_space_removed[mappable_pcode_mask].apply(lambda pcode: postcode_mappings[pcode])
    survey_postcodes.loc[mappable_pcode_idx, 'pcode'] = mapped_postcodes
    survey_postcodes_mapped = postcode_data[postcode_data['Postcode'].isin(survey_postcodes['pcode'])]
    england_postcodes = survey_postcodes_mapped.loc[survey_postcodes_mapped['Country'] == 'England', 'Postcode']
    scotland_postcodes = survey_postcodes_mapped.loc[survey_postcodes_mapped['Country'] == 'Scotland', 'Postcode']
    wales_postcodes = survey_postcodes_mapped.loc[survey_postcodes_mapped['Country'] == 'Wales', 'Postcode']
    NI_postcodes = survey_postcodes_mapped.loc[survey_postcodes_mapped['Country'] == 'Northern Ireland', 'Postcode']

    # scramble england postcodes for confidentiality
    england_postcodes = england_postcodes.reindex(np.random.permutation(england_postcodes.index))

    # send postcodes to a csv for use on the English IMD web api
    england_postcodes.to_csv('data/england_postcodes.csv', index = False, header = False)

    # load in the data generated from the English IMD web api
    england_imd_data = pd.read_excel('data/UK_postcode_IMDs.xlsx', sheet_name = 'english_postcode_IMDs')

    # function for standardising IMD column names
    def rename_imd_cols(df, imd_rank_col, imd_decile_col):
        renamed = df.rename({imd_rank_col: 'IMD rank', imd_decile_col: 'IMD decile'}, axis = 1)
        return renamed

    # filter for IMD columns
    england_imd_columns = ['Index of Multiple Deprivation Rank', 'Index of Multiple Deprivation Decile']
    english_postcode_imds = england_imd_data[['Postcode'] + england_imd_columns]
    english_postcode_imds = rename_imd_cols(english_postcode_imds, *england_imd_columns)

    # load in scotland IMD data and change column names
    scotland_imd_data = pd.read_excel('data/UK_postcode_IMDs.xlsx', sheet_name = 'scottish_postcode_IMDs')
    scotland_imd_columns = ['SIMD2020v2_Rank', 'SIMD2020v2_Decile']
    scottish_postcode_imds = scotland_imd_data.loc[scotland_imd_data['Postcode'].isin(scotland_postcodes), ['Postcode'] + scotland_imd_columns]
    scottish_postcode_imds = rename_imd_cols(scottish_postcode_imds, *scotland_imd_columns)

    # load in wales IMD data and change column names
    wales_imd_data = pd.read_excel('data/UK_postcode_IMDs.xlsx', sheet_name = 'welsh_postcode_IMDs')
    wales_imd_columns = ['WIMD 2019 LSOA Rank', 'WIMD 2019 Overall Decile']
    welsh_postcode_imds = wales_imd_data.loc[wales_imd_data['Welsh Postcode '].isin(wales_postcodes.str.replace(' ', '')), ['Welsh Postcode '] + wales_imd_columns]
    welsh_postcode_imds = rename_imd_cols(welsh_postcode_imds, *wales_imd_columns)
    welsh_postcode_imds.rename({'Welsh Postcode ': 'Postcode'}, axis = 1, inplace = True)
    welsh_postcode_imds.set_index('Postcode', inplace = True)

    # new df for joining postcodes without spacing to postcodes with spaces
    wales_postcode_df = pd.DataFrame(wales_postcodes)
    # set index to be postcode values without spaces
    wales_postcodes_no_spaces = wales_postcodes.str.replace(' ', '')
    wales_postcode_df.set_index(wales_postcodes_no_spaces.values, inplace = True)
    # join on the left index
    wales_postcode_df = wales_postcode_df.merge(welsh_postcode_imds, how = 'left', left_index = True, right_on = 'Postcode')
    # reset index
    welsh_postcode_imds = wales_postcode_df.reset_index(drop = True)

    # loop through northern irish postcodes and get IMD ranks using the get_postcode_rank() function
    pool = Pool(cpu_count())
    postcode_ranks = pool.map(get_postcode_rank, NI_postcodes)
    NI_postcode_imd_ranks = {postcode: rank for postcode, rank in zip(NI_postcodes, postcode_ranks)}
    pool.close()

    # put into data frame
    NI_postcode_imds = pd.DataFrame([(k, v) for k, v in NI_postcode_imd_ranks.items()], columns = ['Postcode', 'IMD rank'])

    # get deciles for each IMD rank
    NI_postcode_imds['IMD decile'] = NI_postcode_imds['IMD rank'].apply(get_decile)

    # concatenate all the imd data
    imd_data_concatenated = pd.concat([english_postcode_imds, scottish_postcode_imds, welsh_postcode_imds, NI_postcode_imds])

    # merge all countries imd data into a single data frame
    survey_imd_data = survey_postcodes.merge(imd_data_concatenated, how = 'left', left_on = 'pcode', right_on = 'Postcode')
    # standardise NA values
    survey_imd_data.loc[survey_imd_data['IMD rank'].isna(), 'IMD rank'] = np.nan
    survey_imd_data.loc[survey_imd_data['IMD decile'].isna(), 'IMD decile'] = np.nan

    # remove extra column
    survey_imd_data.drop('Postcode', axis = 1, inplace = True)
    # save as csv
    survey_imd_data.to_csv('data/Covidence_12Aug20_IMD.csv', index = False)