import pandas as pd
import numpy as np
import re
import math
import os

#######################################################################
#######################################################################

# 1. LOAD OMOP TABLES

# Change directory to where OMOP Vocabulary download is saved
os.chdir('/Users/yc3972/Desktop/DBMI/Courses/GD1_Fall/G4003 Symbolic Methods/Project/off_label/full_OMOP')
# Load all concepts
concept = pd.read_csv('CONCEPT.csv',delimiter='\t',skiprows=1,
                      names = ['concept_id','concept_name',
                              'domain_id','vocabulary_id','concept_class_id','standard_concept',
                              'concept_code','valid_start_date','valid_end_date','invalid_reason'],index_col=False)
# Load all concept relationships
concept_relationship = pd.read_csv('CONCEPT_RELATIONSHIP.csv',delimiter='\t',skiprows=1,
                      names = ['concept_id_1','concept_id_2',
                              'relationship_id',
                               'valid_start_date','valid_end_date','invalid_reason'],
                           index_col=False)

#######################################################################
#######################################################################

# 2. Mapping between vocabs

# 2a. SNOMED Concept_ID --> All NonStandard Concepts
# Get a map of all standard/nonstandard concepts associated with a standrd SNOMED_ID
# Output should match the "Standard to Non-standard map (OMOP)" section on Athena when searching for a given SNOMED_ID 
# Input: single SNOMED Concept ID
# Output: Dataframe with information on the nonStandard concept IDs and their source vocabularies
def get_nonStandard_SNOMED(standard_id):  
    # Get all NonStandard conceptIDs from relationship table identifying "Maps from" for that SNOMED ID 
    NonStandard_Rels = concept_relationship.loc[(concept_relationship['concept_id_1']==standard_id)&
                                                (concept_relationship['relationship_id']=='Mapped from')]  
    # Join nonstandard conceptIDs with OMOP concept table to get data on the concept name, vocabualry source, standard flag, etc. 
    merge_concepts_rels = pd.merge(NonStandard_Rels, concept, 
                                   how="left", left_on='concept_id_2', right_on='concept_id')
    # Pull desired columns, rename columns for clarity, and sort by vocabulary for easy reference 
    concept_map = merge_concepts_rels[['concept_id_1','concept_id_2','concept_name',
                                       'domain_id','vocabulary_id','concept_class_id','standard_concept','concept_code']]
    concept_map =  concept_map.rename(columns = {'concept_id_1':'standard_concept_id','concept_id_':'concept_id'})
    concept_map = concept_map.sort_values(by='vocabulary_id').reset_index(drop=True)
    return concept_map


# 2b. List of ICD Codes --> Associated SNOMED Concept_IDs
# Takes in a df of patient diagnoses from MIMIC data and returns a dictionary of all ICD diagnoses for that patient mapped to standard SNOMED concepts
# Adapted from Annie's Mimic_to_Indication, for full vocabulary set
# Input df with:
    # 'icd_version' 
    # 'icd_code'
# Output: disctionary of each ICD code (including their version) and associated SNOMED Concept_IDs
def diag_ICD_SNOMED(diagnoses):
    std_diag_dict = {}
    for i in range(diagnoses.shape[0]):       
        # Get the ICD code for current diangosis
        version = diagnoses.loc[i, 'icd_version']
        version = 'ICD' + str(version)
        versionCM = str(version) + 'CM'
        code_str = diagnoses.loc[i,'icd_code']
        if len(code_str) < 4:
            code = code_str
        else:
            if (version == 'ICD9') & bool(re.match("E.", code_str)):
                code = code_str[:4] + '.' + code_str[4:]
            #code = '.'.join(code_str[i:i+3] for i in range(0, len(code_str), 3))
            else:
                code = code_str[:3] + '.' + code_str[3:]
        # Get all concept entries for the given ICD code
        icd_match = concept.loc[((concept['vocabulary_id']==version)|(concept['vocabulary_id']==versionCM))&
                                (concept['concept_code']==code)].reset_index()
        # If there is no concept entries for the ICD code, there is no match for that diagnosis in OMOP data
        if icd_match.shape[0] == 0:
            print("No match of " + code + " on OMOP")
        # if both ICD and ICD CM match, then look for exact vocabulary source (i.e. either ICD10 or ICD10CM)
        elif icd_match.shape[0] > 1: 
            icd_concept = icd_match.loc[icd_match.vocabulary_id == version, 'concept_id'].item()

            std_diag = concept_relationship.loc[(concept_relationship['concept_id_1'] == icd_concept) & 
                         (concept_relationship['relationship_id'] == 'Maps to'), 'concept_id_2'].values.tolist()
        # ICD may not match, but ICD CM might match. Use whatever we can
        else:
            icd_concept = icd_match.loc[icd_match.vocabulary_id == versionCM, 'concept_id'].item()

            std_diag = concept_relationship.loc[(concept_relationship['concept_id_1'] == icd_concept) & 
                         (concept_relationship['relationship_id'] == 'Maps to'), 'concept_id_2'].values.tolist()
        given_diag = str(version) + " " + str(code) 
        if given_diag not in std_diag_dict: 
            std_diag_dict[given_diag] = std_diag            
    return std_diag_dict


# 2c. List of ICD Codes --> SNOMED + MeSH Codes
# Takes in a df of patient diagnoses from MIMIC data and returns a dictionary of all ICD diagnoses for that patient mapped to standard SNOMED concepts and NonStandard MeSH Concepts
# Input df with:
    # 'icd_version' 
    # 'icd_code'
    # Uses diag_ICD_SNOMED
# Output: disctionary of each ICD code (including their version) and associated SNOMED + MeSH Concept_IDs
def diag_ICD_MeSH(std_diag_dict):
    
    
    # Initialize ICD -> SNOMED -> MeSH Diagnoses Dictionaries
    mesh_diag_dict = {}
    # Loop through all ICD codes in patient data 
    for entry in std_diag_dict:
        # Intialize a dictionary that will hold the SNOMED standard IDs and Associated MeSH concepts for each ICD code 
        sub_dict = {}
        # Add the SNOMED standard ID
        sub_dict['SNOMED'] = std_diag_dict[entry]
        # Keep track if there is a MeSH concept for that SNOMED ID 
        exists = 0 
        MeSH_list = []
        MeSH_concept_list = []
        # Loop through all SNOMED codes for each ICD code in patient data 
        for SNOMED in std_diag_dict[entry]:
            nonStandard = get_nonStandard_SNOMED(SNOMED) # Get the nonStandard concepts for that 
            if(nonStandard.loc[(nonStandard['vocabulary_id']=='MeSH') & (nonStandard['concept_class_id']=='Main Heading')].empty)==False:
                exists = 1
                MeSH_list.extend(nonStandard.loc[(nonStandard['vocabulary_id']=='MeSH') & 
                                                 (nonStandard['concept_class_id']=='Main Heading')]['concept_id_2'].values.tolist())
                MeSH_concept_list.extend(nonStandard.loc[(nonStandard['vocabulary_id']=='MeSH') & 
                                                         (nonStandard['concept_class_id']=='Main Heading')]
                                         ['concept_name'].values.tolist())        
        sub_dict['MeSH'] = MeSH_list
        sub_dict['MeSH_concept'] = MeSH_concept_list
        if(exists): 
            mesh_diag_dict[entry] = sub_dict              
    return mesh_diag_dict


# 2d. Single MeSH term --> Any other Vocab
# Get the Concept_ID from OMOP associated with the MeSH Term
def get_MeSH_conceptID(MeSH_term):
    MeSH_concept = concept.loc[(concept['vocabulary_id']=='MeSH')&(concept['concept_name']==MeSH_term)]
    if(len(MeSH_concept)>0):
        MeSH_conceptId = MeSH_concept.reset_index()['concept_id'][0]
        return MeSH_conceptId
    else:
        return

# Get a map of all standard/nonstandard concepts associated with that MeSH term
def get_nonStandard_MeSH(MeSH_term):
    # Identify concept_id for MeSH term 
    MeSH_conceptId = get_MeSH_conceptID(MeSH_term)
    
    # Get the entry in the CONCEPT_RELATIONSHIP table for that MeSH term
    MeSH_rel_entry = concept_relationship.loc[concept_relationship['concept_id_1']==MeSH_conceptId].reset_index() # NonStandard --> Standard  
    if(len(MeSH_rel_entry)>0):
        # Get the SNOMED standard ID for that concept
        standard_id = MeSH_rel_entry['concept_id_2'][0]
    
        # Get all NonStandard conceptIDs from relationship tables for that SNOMED ID 
        NonStandard_Concepts = concept_relationship.loc[(concept_relationship['concept_id_1']==standard_id)&(concept_relationship['relationship_id']=='Mapped from')]  # Get all NonStandard concept IDs


        # Join nonstandard conceptIDs with OMOP concept table to get data on the concept name, vocabualry source, standard/non, etc. 
        merge_concepts_rels = pd.merge(NonStandard_Concepts, concept,
                                how="left", left_on='concept_id_2', right_on='concept_id')

        # Pull desired columns, rename columns for clarity, and sort by vocabulary for easy reference 
        concept_map = merge_concepts_rels[['concept_id_1','concept_id_2','concept_name','domain_id','vocabulary_id','concept_class_id','standard_concept','concept_code']]
        concept_map =  concept_map.rename(columns = {'concept_id_1':'standard_concept_id','concept_id_':'concept_id'})
        concept_map = concept_map.sort_values(by='vocabulary_id').reset_index(drop=True)

        return concept_map
    else:
        return


# Get the concept codes associated with a MeSH term from a specified alternative vocabulary (i.e. ICD10/ICD10CM) 
def get_MeSH_vocab_codes(MeSH_term,vocab):
    concept_map = get_nonStandard_MeSH(MeSH_term)
    if isinstance(concept_map, pd.DataFrame):
        codes = (concept_map.loc[concept_map['vocabulary_id']==vocab])['concept_code'].values.tolist()
        return codes
    else:
        return



#######################################################################
#######################################################################

# 3. Parallel + Hierarchical SNOMED Traversals
 
# 3a. Standard SNOMED -> NonStandard SNOMED
### Get All related SNOMED codes for given SNOMED code, HORIZONTALLY
### For mapping, refer to vocab_mapping.py, get_nonStandard_SNOMED(j)
### We specifically want to get all the SNOMED codes that map to out standard, for efficient comparison with indication
def get_related(std_dict):
    related = {}
    for icd in std_dict.keys():
        std_codes = std_dict[icd]
        for j in std_codes:
            all_related = get_nonStandard_SNOMED(j)
            all_related = all_related.loc[all_related.vocabulary_id == 'SNOMED','concept_id_2'].values.tolist()
            if icd not in related:
                related[icd] = all_related
            else:
                related[icd].extend(all_related)
    return related


# 3b. Standard SNOMED --> All Parent/Child Relationships
def get_hierarchy(std_diag_dict):    
    snomed_hierarchy = {}
    snomed_parent = {}
    snomed_child = {}
    # Loop through all ICD codes in patient data 
    for entry in std_diag_dict:
        # Intialize a dictionary that will hold the SNOMED standard IDs and Associated MeSH concepts for each ICD code 
        sub_dict = {}
        # Add the SNOMED standard ID
        sub_dict['SNOMED'] = std_diag_dict[entry]
        # Loop through all SNOMED codes for each ICD code in patient data 
        parents = []
        children = []
        for SNOMED in std_diag_dict[entry]:
            parents_temp = concept_relationship.loc[(concept_relationship['concept_id_1']==SNOMED)&
                                                (concept_relationship['relationship_id']=='Is a')]['concept_id_2'].values.tolist()
            parents = parents + parents_temp
            children_temp = concept_relationship.loc[(concept_relationship['concept_id_1']==SNOMED)&
                                                (concept_relationship['relationship_id']=='Subsumes')]['concept_id_2'].values.tolist()
            children = children + children_temp
        sub_dict['SNOMED_parents'] = parents
        sub_dict['SNOMED_children'] = children
        snomed_hierarchy[entry] = sub_dict 
        snomed_parent[entry] = std_diag_dict[entry] + parents
        snomed_child[entry] = std_diag_dict[entry] + children
    return snomed_hierarchy,snomed_parent,snomed_child

#######################################################################
#######################################################################

# 4. RxNorm, NDF-RT Mappings: extract and map drug concepts from non-standard terms and vice versa. Used to convert from NDC to RxNorm ingredient, and to NDF-RT codes.

# Drug dictionary that stores the standard RxNorm codes for the given drugs
# Used if we were to traverse across drugs (to find related other ndc codes)
# NDC codes --> RxNorm standard 
def NDC_RxNorm(drugslist):
    std_drug_dict = {}
    for i in range(drugslist.shape[0]):
        drug = drugslist.loc[i,'drug']
        ndc = drugslist.loc[i, 'ndc']

        # ndc = 0 means there is no ndc code
        if ndc == '0':
            continue
        ndc_id = concept.loc[concept.concept_code == ndc, 'concept_id'].values.tolist()

        if len(ndc_id) == 0: # there is no matching standard
            continue
        ndc_id2 = ndc_id[0]
        standard_id = maps_to(ndc_id2)
        ingredient = find_ingredient(drug)
        # if rxnorm ingredient cannot be found directly, manually find via relationships
        if len(ingredient) == 0 or len(ingredient) > 1: 
            ingredient = find_ingredient2(standard_id)

        if drug not in std_drug_dict: 
            std_drug_dict[drug] = ingredient
    return std_drug_dict

# Find RxNorm ingredient of given drug based on OMOP standard ID, based on relationship type
# If you had one dataset for both drug and diagnosis, change "drug_mapping" to the name of the relationship set
def find_ingredient2(standard_id):
    relationships = concept_relationship.loc[(concept_relationship.concept_id_1 == standard_id),'relationship_id'].unique()
    
    if 'RxNorm has ing' in relationships: 
        ingredient = concept_relationship.loc[(concept_relationship.concept_id_1 == standard_id) & 
                                  (concept_relationship.relationship_id == 'RxNorm has ing'), 'concept_id_2'].values.tolist()
        return ingredient
    
    elif 'Has precise ing' in relationships:
        ingredient = concept_relationship.loc[(concept_relationship.concept_id_1 == standard_id) & 
                                  (concept_relationship.relationship_id == 'Has precise ing'), 'concept_id_2'].values.tolist()
        return ingredient
    
    elif 'Has precise ingredient' in relationships:
        ingredient = concept_relationship.loc[(concept_relationship.concept_id_1 == standard_id) & 
                                  (concept_relationship.relationship_id == 'Has precise ingredient'), 'concept_id_2'].values.tolist()
        return ingredient
    
    elif 'Has brand name' in relationships:
        ingredient = concept_relationship.loc[(concept_relationship.concept_id_1 == standard_id) & 
                                  (concept_relationship.relationship_id == 'Has brand name'), 'concept_id_2'].values.tolist()
        ingredient2 = ingredient[0]
        return find_ingredient2(ingredient2)
    
    elif 'Consists of' in relationships:
        composition = concept_relationship.loc[(concept_relationship.concept_id_1 == standard_id) & 
                                  (concept_relationship.relationship_id == 'Consists of'), 'concept_id_2'].values
        if len(composition) > 1:
            print("more than one consists of relationships")
        else:
            composition = composition.item()
            return find_ingredient2(composition)
    else: 
        ingredient = standard_id
        return [ingredient] 

# Find RxNorm ingredient of drug based on concept class id
def find_ingredient(drug):
    all_ingredients = concept.loc[(concept.concept_class_id == 'Ingredient') & (concept.vocabulary_id == 'RxNorm')]
    ingredient = all_ingredients.loc[all_ingredients.concept_name.str.contains(drug, case = False, na = False, regex = False), 'concept_id'].values.tolist()
    return ingredient


 # Traverse drug mapping until there is no more maps to
 # feed in either drug mapping or diagnosis mapping data (if both are combined, just feed in relationship data)
def maps_to(conceptid):
    mapping = concept_relationship
    relationships = mapping.loc[mapping.concept_id_1 == conceptid, 'relationship_id'].unique()
    if 'Maps to' not in relationships:
        return conceptid
    
    else:
        mapped = mapping.loc[(mapping.concept_id_1 == conceptid) & 
                                        (mapping.relationship_id == 'Maps to'), 'concept_id_2'].values.tolist()
        
        root_concept = mapped[0]
        if root_concept != conceptid:
            newconcept = root_concept
            maps_to(newconcept)
        
        return root_concept

# Find all relevant (non-standard) concepts that map to the standard concept
# feed in either drug mapping or diagnosis mapping data (if both are combined, just feed in relationship data)
def mapped_from(conceptid):
    mapping = concept_relationship
    relationships = mapping.loc[mapping.concept_id_1 == conceptid, 'relationship_id'].unique()
    if 'Mapped from' not in relationships:
        return conceptid
    else:
        all_codes = mapping.loc[(mapping.concept_id_1 == conceptid) & 
                                   (mapping.relationship_id == 'Mapped from')]
        return all_codes



# Find NDFRT terms equivalent 
# Takes in rxnorm ingredient ID, returns NDFRT equivalent of that concept
def find_NDFRT(rxnorm_id):
    #check if there is a rxnorm-ndfrt mapping
    relationships = concept_relationship.loc[concept_relationship.concept_id_1 == rxnorm_id, 'relationship_id'].unique()
    if 'RxNorm - NDFRT eq' not in relationships:
        pass
    else: 
        ndfrt = concept_relationship.loc[(concept_relationship.concept_id_1 == rxnorm_id) & 
                                 (concept_relationship.relationship_id == 'RxNorm - NDFRT eq'), 'concept_id_2'].values.tolist()
        return ndfrt

    
# Finding the NDFRT equivalent of that rxnorm ingredient, necessary for the treat relationships
# Find mapping code in find_NDFRT
# Takes in rxnorm ingredient list, returns dictionary of ndfrt equivalents
def all_ndfrt(rxnorm_list):
    ndfrt_list = {}
    for drug in rxnorm_list.keys():
        if rxnorm_list[drug] is None:
            print(drug, rxnorm_list[drug])
            continue
        rxnorm = rxnorm_list[drug][0]
        ndfrt = find_NDFRT(rxnorm)
        if ndfrt == None:
            print(drug, "doesn't have ndfrt equivalent")
            continue
        if drug not in ndfrt_list:
            ndfrt_list[drug] = ndfrt
    return ndfrt_list
    
# Utilize the same nonstandard mapping, but we want NDFRT codes instead of SNOMED
# Takes in standard dictionary of drugs, returns all horizontal nonstandard mapping of NDFRT values
def all_related_NDFRT(std_dict):
    related = {}
    for icd in std_dict.keys():
        std_codes = std_dict[icd]
        for j in std_codes:
            all_related = get_nonStandard_SNOMED(j)
            all_related = all_related.loc[all_related.vocabulary_id == 'NDFRT','concept_id_2'].values.tolist()
            if icd not in related:
                related[icd] = all_related
            else:
                related[icd].extend(all_related)
    return related

# Check treat by relationships in OMOP using NDF_RT data
# Takes in a list of NDFRT drug codes and returns all the NDFRT diagnosis codes that may be treated by the drug
def omop_indications(ndfrt_list):
    indications = []
    for i in range(len(ndfrt_list)):
        concepts = concept_relationship.loc[(concept_relationship.concept_id_1 == ndfrt_list[i]) & 
                         (concept_relationship.relationship_id == 'May treat'), 'concept_id_2'].values.tolist()
        indications.extend(concepts)
    return indications

# Now, we try to map these NDFRT indications back to SNOMED so that we can compare with our diagnoses
# Takes in list of ndfrt indications list and returns snomed list
def ndfrt_SNOMED(ndfrt_list):
    snomed_list = []
    for ndfrt in ndfrt_list:
        snomed = concept_relationship.loc[(concept_relationship.concept_id_1 == ndfrt)&
                        (concept_relationship.relationship_id == 'Ind/CI - SNOMED'), 'concept_id_2'].values.tolist()
        snomed_list.extend(snomed)
    return snomed_list
