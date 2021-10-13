import re
import xml.etree.ElementTree as ET

tree = ET.parse('../data/drugbank.xml')
root = tree.getroot()

ns = '{http://www.drugbank.ca}'

# pattern for newline characters in drug names
newline_pattern = re.compile('^\s*\\n\s*$')
# pattern for parentheses in drug names (suffixes)
parentheses_pattern = re.compile('\(+.*\)+\s*$')
# pattern for dosages included in drug names
dosage_pattern = re.compile('[\d.%]+\s*\w+/*\d*[\w.\s]*\s*$')
# pattern for non alphanumeric characters at the end or beginning of string
punctuation_pattern = re.compile('^[^\w]+|[^\w]+$')
# dictionary to store all drug names with pointer to the drugbank id
drug_dictionary = {}
db_id_dictionary = {}
# list for storing canonical drug names - these should not be overwritten in the dictionary
canonical_names = []

# loop through drug entries
for drug in root:

    # get entry name and drugbank ID
    drug_id = drug.findtext(ns + "drugbank-id[@primary='true']")
    name = drug.findtext(ns + 'name').lower()
    canonical_names.append(name)
    drug_dictionary[name] = set([drug_id])

    # add drug aliases
    international_brands = {elem.text.lower() for elem in drug.findall('{ns}international-brands/{ns}international-brand/{ns}name'.format(ns = ns))}
    synonyms = {elem.text.lower() for elem in drug.findall('{ns}synonyms/{ns}synonym'.format(ns=ns))}
    products = {elem.text.lower() for elem in drug.findall('{ns}products/{ns}product/{ns}name'.format(ns = ns))}
    aliases = international_brands.union(synonyms, products)

    # trim suffix (an ending phrase contained in parentheses)
    aliases_suffix_trimmed = {parentheses_pattern.sub('', entry) for entry in aliases if not newline_pattern.search(entry)}
    # trim dosage
    aliases_dosage_removed = {dosage_pattern.sub('', entry) for entry in aliases_suffix_trimmed}
    # trim punctuation
    alias_punct_removed = {punctuation_pattern.sub('', entry) for entry in aliases_dosage_removed}
    # remove empty entries
    aliases_cleaned = {entry for entry in alias_punct_removed if entry}

    # add to the dictionary
    for alias in aliases_cleaned:
        # if an alias for another drug is already in the canonical names, don't change and just continue
        if alias in canonical_names:
            continue
        # otherwise if its an alias in the drug dictionary, take the union with the existing entry
        elif alias in drug_dictionary:
            drug_dictionary[alias] = drug_dictionary[alias].union({drug_id})
        # otherwise make a new entry
        else:
            drug_dictionary[alias] = set([drug_id])
        # add the id to the id dictionary, paired to alias
        if drug_id in db_id_dictionary:
            db_id_dictionary[drug_id] = db_id_dictionary[drug_id].union({alias})
        else:
            db_id_dictionary[drug_id] = set([alias])


# add drug mixture products
mixture_dict = {}      
for drug in root:

    # get mixture names and ingredients from the xml tree
    mixture_names = drug.findall('{ns}mixtures/{ns}mixture/{ns}name'.format(ns=ns))
    mixture_ingredients = drug.findall('{ns}mixtures/{ns}mixture/{ns}ingredients'.format(ns=ns))

    # loop through names and indgredients
    for name, ingredients in zip(mixture_names, mixture_ingredients):

        # only map mixture products with 2 or more ingredients
        if len(ingredients.text.split('+')) > 1:
            # filter some common patterns
            name_suffix_trimmed = parentheses_pattern.sub('', name.text.lower())
            name_cleaned = dosage_pattern.sub('', name_suffix_trimmed).strip()

            # if name is already in the canonical names list, move to the next mixture
            if name_cleaned in canonical_names:
                continue
            # get set of ingredients
            ingredients_set = {ingredient.lower().strip() for ingredient in ingredients.text.split('+')}
            mapped_db_ids = set()

            for ingredient in ingredients_set:
                # if the ingredient is in the dictionary, add to the mapped db ids
                if ingredient in drug_dictionary:
                    mapped_db_ids = mapped_db_ids.union(drug_dictionary.get(ingredient))

            # if the name is already in the mixture dictionary, take the union
            if name_cleaned in mixture_dict:
                mixture_dict[name_cleaned] = mixture_dict[name_cleaned].union(mapped_db_ids)

            # otherwise, if there are mapped drugbank ids, add to the mixture dictionary under the name
            elif mapped_db_ids:
                mixture_dict[name_cleaned] = mapped_db_ids
                
# go through the mixture dictionary and add entries to the drug dictionary
for mixture in mixture_dict:
    if mixture in drug_dictionary:
        drug_dictionary[mixture] = drug_dictionary[mixture].union(mixture_dict[mixture])
    else:
        drug_dictionary[mixture] = mixture_dict[mixture]

## adding first word names to the drug dictionary ##

# list of first names that have been added
added_first_names = []

# check the drug dictionary to see if the first word of each entry is a separate entry
# if not save the first word of the name to the drug dictionary mapped to the drugbank ids of the full name
for drug_name in list(drug_dictionary):

    name_split = drug_name.split(' ')

    if len(name_split) > 1:
        first_word = name_split[0]
        # if the first word is in the drug dictionary, check that it has not been added in this loop
        if first_word in drug_dictionary:
            # check for ambiguity - i.e. if the first name is already added and is different to another potential mapping
            if first_word in added_first_names and drug_dictionary[first_word] != drug_dictionary[drug_name]:
                drug_dictionary.pop(first_word)
            else:
                continue

        # if the first word is not already in the drug dictionary, add it
        else:
            drug_dictionary[first_word] = drug_dictionary[drug_name]
            # add to list to track names that have been added
            added_first_names.append(first_word)