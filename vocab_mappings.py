import pandas as pd
import numpy as np
import re
import math
import os

#######################################################################
#######################################################################

# 1. LOAD OMOP TABLES

# Change directory to where OMOP Vocabulary download is saved
os.chdir('/Users/khs2138/OneDrive - cumc.columbia.edu/Symbolic_Project/OMOP_vocabulary')
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
 
# To be modified for compatibility with OMOP full dataset:

# MIMIC stores drug data as NDC codes
drug_concept = pd.read_csv('./data/DRUG_VOCAB/CONCEPT.csv', sep = '\t', header = 0)
# What is the corresponding RxNorm standard code for that drug?
drug_mapping = pd.read_csv('./data/DRUG_VOCAB/CONCEPT_RELATIONSHIP.csv', sep = '\t', header = 0)
# ICD_SNOMED mapping
diag_concept = pd.read_csv('./data/SNOMED_ICD/CONCEPT.csv', sep = '\t', header = 0)
diag_mapping = pd.read_csv('./data/SNOMED_ICD/CONCEPT_RELATIONSHIP.csv', sep = '\t', header = 0)

# Find RxNorm ingredient of given drug based on OMOP standard ID, based on relationship type
# If you had one dataset for both drug and diagnosis, change "drug_mapping" to the name of the relationship set
def find_ingredient2(standard_id):
    relationships = drug_mapping.loc[(drug_mapping.concept_id_1 == standard_id),'relationship_id'].unique()
    
    if 'RxNorm has ing' in relationships: 
        ingredient = drug_mapping.loc[(drug_mapping.concept_id_1 == standard_id) & 
                                  (drug_mapping.relationship_id == 'RxNorm has ing'), 'concept_id_2'].item()
        return ingredient
    
    elif 'Has precise ing' in relationships:
        ingredient = drug_mapping.loc[(drug_mapping.concept_id_1 == standard_id) & 
                                  (drug_mapping.relationship_id == 'Has precise ing'), 'concept_id_2'].item()
        return ingredient
    
    elif 'Has precise ingredient' in relationships:
        ingredient = drug_mapping.loc[(drug_mapping.concept_id_1 == standard_id) & 
                                  (drug_mapping.relationship_id == 'Has precise ingredient'), 'concept_id_2'].item()
        return ingredient
    
    elif 'Has brand name' in relationships:
        ingredient = drug_mapping.loc[(drug_mapping.concept_id_1 == standard_id) & 
                                  (drug_mapping.relationship_id == 'Has brand name'), 'concept_id_2'].item()
        return find_ingredient2(ingredient)
    
    elif 'Consists of' in relationships:
        composition = drug_mapping.loc[(drug_mapping.concept_id_1 == standard_id) & 
                                  (drug_mapping.relationship_id == 'Consists of'), 'concept_id_2'].values
        if len(composition) > 1:
            print("more than one consists of relationships")
        else:
            composition = composition.item()
            return find_ingredient2(composition)
    else: 
        ingredient = standard_id
        return ingredient 

# Find RxNorm ingredient of drug based on concept class id
def find_ingredient(drug):
    all_ingredients = drug_concept.loc[drug_concept.concept_class_id == 'Ingredient']
    ingredient = all_ingredients.loc[all_ingredients.concept_name.str.contains(drug, case = False, na = False, regex = False), 'concept_id'].values.tolist()
    return ingredient


 # Traverse drug mapping until there is no more maps to
 # feed in either drug mapping or diagnosis mapping data (if both are combined, just feed in relationship data)
def maps_to(mapping_data, conceptid):
    mapping = mapping_data
    relationships = mapping.loc[mapping.concept_id_1 == conceptid, 'relationship_id'].unique()
    if 'Maps to' not in relationships:
        return conceptid
    
    else:
        root_concept = mapping.loc[(mapping.concept_id_1 == conceptid) & 
                                        (mapping.relationship_id == 'Maps to'), 'concept_id_2'].item()
        if root_concept != conceptid:
            newconcept = root_concept
            maps_to(mapping, newconcept)
        
        return root_concept

# Find all relevant (non-standard) concepts that map to the standard concept
# feed in either drug mapping or diagnosis mapping data (if both are combined, just feed in relationship data)
def mapped_from(mapping_data, conceptid):
    mapping = mapping_data
    relationships = mapping.loc[mapping.concept_id_1 == conceptid, 'relationship_id'].unique()
    if 'Mapped from' not in relationships:
        return conceptid
    else:
        all_codes = mapping.loc[(mapping.concept_id_1 == conceptid) & 
                                   (mapping.relationship_id == 'Mapped from')]
        return all_codes



# Find NDFRT terms equivalent 
def find_NDFRT(rxnorm_id):
    #check if there is a rxnorm-ndfrt mapping
    relationships = drug_mapping.loc[drug_mapping.concept_id_1 == rxnorm_id, 'relationship_id'].unique()
    if 'RxNorm - NDFRT eq' not in relationships:
        pass
    else: 
        ndfrt = drug_mapping.loc[(drug_mapping.concept_id_1 == rxnorm_id) & 
                                 (drug_mapping.relationship_id == 'RxNorm - NDFRT eq'), 'concept_id_2'].values.tolist()
        return ndfrt

