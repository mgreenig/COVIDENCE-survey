from urllib.request import Request, urlopen
from bs4 import BeautifulSoup as bs
import re
from string import ascii_uppercase
from math import ceil
import pickle

# get the DrugBank drug dictionary
from Parse_drugbank import drug_dictionary

print('DrugBank XML tree parsed, pulling compounds from EMC...')

## functions for parsing html data from EMC ##

# get the html data from a url
def get_html(url):
    req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    webpage = urlopen(req).read()
    html = bs(webpage, 'html.parser')
    return html

# regex pattern for the number of drug results under a certain letter
n_results_pattern = re.compile('\d+(?= results found)')
# function for getting the number of drug results under a certain letter
def get_n_results(html):
    n_results_text = html.find('span', {'class': 'search-paging-view'}).text
    n_results = int(n_results_pattern.search(n_results_text).group(0))
    return n_results

# get drugs and corresponding links
def get_links_on_page(url):
    html = get_html(url)
    drugs = [elem.find('h2').text.strip('\n').lower() for elem in html.find_all('div', {'class': 'row data-row'})]
    links = [elem.find('h2').find('a').get('href') for elem in html.find_all('div', {'class': 'row data-row'})]
    # only save keys corresponding to the lowercase first word of every entry
    link_dict = {drug: link for drug, link in zip(drugs, links)}
    # only save drugs that are not found in the drug dictionary
    link_dict_filtered = {drug: link for drug, link in link_dict.items() if
                          re.sub('[^(\w|/|\-)].*$', '', drug).lower() not in drug_dictionary}
    return link_dict_filtered


# get the active ingredients from the link to a drug on EMC
def get_active_ingredients(link):
    drug_url = 'https://www.medicines.org.uk' + link
    try:
        drug_html = get_html(drug_url)
        active_ingredients = [elem.text.strip().lower() for elem in drug_html.find('div', {'class': 'col-xs-12 col-sm-6'}).find_all('li')]
        # split individual ingredient listings containing semi-colons
        active_ingredients = [ingr for ingrs in active_ingredients for ingr in ingrs.split('; ')]
    # return none if an error is thrown
    except:
        print('Error on {}'.format(link))
        active_ingredients = None
    return active_ingredients

if __name__ == '__main__':

    # list for all urls we need to pull active ingredients from
    all_urls = []

    # loop through letters and get the drug names under each letter
    for letter in ascii_uppercase:
        # get the number of pages under each letter
        url = 'https://www.medicines.org.uk/emc/browse-medicines/{}'.format(letter)
        html = get_html(url)
        n_results = get_n_results(html)
        n_pages = ceil(n_results/200)
        # get the all the drug pages for the letter and append to list
        letter_urls = ['https://www.medicines.org.uk/emc/browse-medicines?prefix={}&offset={}&limit=200'.format(letter, (i*200)+1) for i in range(n_pages)]
        all_urls.extend(letter_urls)


    drug_links = map(get_links_on_page, all_urls)

    # dictionary for saving drug links
    all_drug_links = {}

    # add the links to the global dictionary
    for link_dict in drug_links:
        all_drug_links.update(link_dict)

    print('Links for unmapped compounds added, getting active ingredients...')

    # dictionary for storing active ingredients
    all_active_ingredients = {}

    # loop through drug links and get active ingredients
    for drug, link in all_drug_links.items():
        active_ingredients = get_active_ingredients(link)
        # for active ingredients not in the drug dictionary, take the first word only
        active_ingredients = [ingr if ingr in drug_dictionary else re.sub('\s+.*$', '', ingr) for ingr in active_ingredients]
        all_active_ingredients[drug] = active_ingredients

    # dictionary drug links where get_active_ingredients() returned None
    unmapped = {drug: link for drug, link in all_drug_links.items() if all_active_ingredients[drug] is None}

    # map the unmapped drugs
    for drug, link in unmapped.items():
        all_active_ingredients[drug] = get_active_ingredients(link)

    # filter drugs that returned None
    all_active_ingredients = {drug: ingrs for drug, ingrs in all_active_ingredients.items() if ingrs}

    print('Active ingredients pulled, adding to drug dictionary...')

    # dictionary for drugbank IDs mapped to drugs
    drug_ingredient_IDs = {}

    # list for drug names with ambiguous active ingredients
    ambiguous_drug_names = []

    # loop through drugs in the dictionary and save active ingredient IDs under shortened names
    for drug, ingredients in all_active_ingredients.items():
        # get all characters before first non-alphanumeric character including / and -
        drug_shortened = re.sub('[^(\w|/|\-)]+.+$', '', drug)
        # get drug dictionary IDs for ingredients
        active_ingredient_IDs = set().union(*[drug_dictionary.get(ingr) for ingr in ingredients if drug_dictionary.get(ingr)])
        # if the drug is already in the dictionary and the existing value does not match the active ingredient IDs, add to ambiguous list
        if drug_shortened in drug_ingredient_IDs and drug_ingredient_IDs[drug_shortened] != active_ingredient_IDs:
            ambiguous_drug_names.append(drug_shortened)
        # otherwise if the drug is not present in the dictionary, add it
        elif drug_shortened not in drug_ingredient_IDs:
            drug_ingredient_IDs[drug_shortened] = active_ingredient_IDs

    # remove ambiguous entries from the dictionary
    drug_ingredient_IDs = {drug: ingrs for drug, ingrs in drug_ingredient_IDs.items() if drug not in ambiguous_drug_names}

    # loop through active ingredient mappings
    for drug in drug_ingredient_IDs:
        # if drugbank ids were obtained for the drug's active ingredients, add to the drug dictionary
        if drug_ingredient_IDs[drug]:
            drug_dictionary[drug] = drug_ingredient_IDs[drug]

    # save the dictionary to a pickle
    pickle.dump(drug_dictionary, open('data/drug_dictionary.p', 'wb'))