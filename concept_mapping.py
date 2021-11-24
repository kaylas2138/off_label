import pandas as pd
import os


# get your current working directory
os.getcwd()
# set your working directory to the directory in which you have all your mimic data files
os.chdir('/Users/yc3972/Desktop/DBMI/Courses/GD1_Fall/G4003 Symbolic Methods/Project/off_label')

"""Reading in datasets:
I've downloaded OMOP vocabulary separately for drugs and diagnoses; if you have one vocabulary for all,
modify the code to read only the concept and relationships files"""

# MIMIC stores drug data as NDC codes
drug_concept = pd.read_csv('./data/DRUG_VOCAB/CONCEPT.csv', sep = '\t', header = 0)
# What is the corresponding RxNorm standard code for that drug?
drug_mapping = pd.read_csv('./data/DRUG_VOCAB/CONCEPT_RELATIONSHIP.csv', sep = '\t', header = 0)
# ICD_SNOMED mapping
diag_concept = pd.read_csv('./data/SNOMED_ICD/CONCEPT.csv', sep = '\t', header = 0)
diag_mapping = pd.read_csv('./data/SNOMED_ICD/CONCEPT_RELATIONSHIP.csv', sep = '\t', header = 0)

# As an example for Kayla:
# %cd vocabulary_download_v5_{a57bccfe-1f01-4e58-aed1-3284504a7f79}_1636826231531/
# concept = pd.read_csv('CONCEPT.csv',delimiter='\t',skiprows=1,
#                       names = ['concept_id','concept_name',
#                               'domain_id','vocabulary_id','concept_class_id','standard_concept',
#                               'concept_code','valid_start_date','valid_end_date','invalid_reason'],index_col=False)
# concept_relationship = pd.read_csv('CONCEPT_RELATIONSHIP.csv',delimiter='\t',skiprows=1,
#                       names = ['concept_id_1','concept_id_2',
#                               'relationship_id',
#                                'valid_start_date','valid_end_date','invalid_reason'],
#                            index_col=False)


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
                                   (mapping.relationship_id == 'Mapped from'), 'concept_id_2'].values.tolist()
        return all_codes



# Find NDFRT terms equivalent 
def find_NDFRT(rxnorm_id):
    #check if there is a rxnorm-ndfrt mapping
    relationships = drug_mapping.loc[drug_mapping.concept_id_1 == rxnorm_id, 'relationship_id'].unique()
    if 'RxNorm - NDFRT eq' not in relationships:
        print("there is no NDF_RT equivalent to this RxNorm concept")
    else: 
        ndfrt = drug_mapping.loc[(drug_mapping.concept_id_1 == rxnorm_id) & 
                                 (drug_mapping.relationship_id == 'RxNorm - NDFRT eq'), 'concept_id_2'].values.tolist()
        return ndfrt