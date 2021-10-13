from urllib.request import urlopen
from urllib.parse import quote
from bs4 import BeautifulSoup as bs
import pandas as pd
import re

# if running the file, pull new data
if __name__ == '__main__':

    bnf_drug_url = urlopen('https://bnf.nice.org.uk/drug/')

    bnf_drugs_html = bs(bnf_drug_url, 'html.parser')

    # drug names are kept as a 'span' html object
    drug_html_data = bnf_drugs_html.find_all('span')

    # clean the text
    drug_names = [elem.text.lower().strip() for elem in drug_html_data if elem.text]
    # remove last element, which does not correspond to a drug
    drug_names.pop()

    # dictionary for storing the BNF classification for each drug
    drug_classifications = {}

    # loop through drugs, access url, and get the classification for each
    for drug in drug_names:
        drug_name_for_url = re.sub('[\s,]+', '-', drug)
        drug_name_for_url = re.sub('[^\w-]', '', drug_name_for_url)
        drug_name_for_url = quote(drug_name_for_url)
        drug_url = urlopen('https://bnf.nice.org.uk/drug/{}.html'.format(drug_name_for_url))
        drug_html = bs(drug_url, 'html.parser')
        primary_classification_html_obj = drug_html.find_all('a', {'class': 'classification primary-classification'})
        secondary_classification_html_obj = drug_html.find_all('a', {'class': 'classification secondary-classification'})
        # return the list of primary classifications if it exists, else 'None'
        if len(primary_classification_html_obj) > 0:
            primary_classification = [elem.text.strip() for elem in primary_classification_html_obj]
        # save None as a list so we can use str.join later
        else:
            primary_classification = ['None']
        # same with secondary classification
        if len(secondary_classification_html_obj) > 0:
            secondary_classification = [elem.text.strip() for elem in secondary_classification_html_obj]
        else:
            secondary_classification = ['None']
        drug_classifications[drug] = [primary_classification, secondary_classification]

    bnf_classes = pd.DataFrame(drug_classifications).T
    bnf_classes.columns = ['primary', 'secondary']

    # split multi-drug formulations into multiple drugs
    bnf_classes = bnf_classes.rename_axis('drugs').reset_index()
    bnf_classes['drugs'] = bnf_classes['drugs'].str.split('\s+with\s+|\s+and\s+')
    bnf_classes['drugs'] = bnf_classes['drugs'].apply(lambda drugs: [drug.strip() for drug in drugs])

    # join columns that contain lists/sets
    bnf_classes['drugs'] = bnf_classes['drugs'].str.join('; ')
    bnf_classes['primary'] = bnf_classes['primary'].str.join('; ')
    bnf_classes['secondary'] = bnf_classes['secondary'].str.join('; ')

    bnf_classes.to_csv('../data/bnf_drug_classifications.csv', index_label='entry')