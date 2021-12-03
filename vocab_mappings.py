
import pandas as pd
import numpy as np
import re
import math
import os

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
def diag_ICD_MeSH(diagnoses):
    # Map ICD --> SNOMED
    std_diag_dict = diag_ICD_SNOMED(diagnoses)     
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
            if((nonStandard.loc[nonStandard['vocabulary_id']=='MeSH'].empty) == False):
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


# In[15]:


# Get the concept codes associated with a MeSH term from a specified alternative vocabulary (i.e. ICD10/ICD10CM) 
def get_MeSH_vocab_codes(MeSH_term,vocab):
    concept_map = get_nonStandard_MeSH(MeSH_term)
    if isinstance(concept_map, pd.DataFrame):
        codes = (concept_map.loc[concept_map['vocabulary_id']==vocab])['concept_code'].values.tolist()
        return codes
    else:
        return





