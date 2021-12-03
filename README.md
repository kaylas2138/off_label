# off_label
Automated mapping of Electronic Health Records to document off-label drug use 

off_label.py >> Extracts relevant data from MIMICIV publicly available clinical data, which can be accessed at https://physionet.org/content/mimiciv/0.4/. 
vocab_mappings.py >> Python script that defines functions that extract standard to non-standard terms or vice versa. Used to convert from ICD to SNOMED, and to MeSH codes.
concept_mapping.py >> Python script used to extract and map drug concepts from non-standard terms and vice versa. Used to convert from NDC to RxNorm ingredient, and to NDF-RT codes. 

**We start by extracting formal indication dat from Drug Central, in our script DrugInd_Extraction.ipynb. We also leverage Clinical Trial data from AACT. Please refer to link within script for more details on each dataset.
We extract clinical data from MIMIC, and examine the list of diagnoses and drugs given to a patient during a specific visit. We map the nonstandard codes for both diagnoses and drugs to their standard equivalent (SNOMED and RxNorm ingredient, respectively). 
We compare our diagnoses list with known indications from the drugs based on three indications sources: Drug Central, AACT, and NDF-RT. 
 
We analyze implications of our findings in Implications.ipynb. 
