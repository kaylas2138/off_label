# Project Description

Automated mapping of Electronic Health Records to document off-label drug use. There are thousands of off-label medications that are being prescribed to patients on a daily basis. For this project, we have developed an automated method to detect and analyze off-label prescribing patterns by comparing diagnosis codes to FDA-approved indications, commonly occurring off-label uses, and uses being studied in clinical trials for prescribed medications within the [MIMICIV dataset](https://physionet.org/content/mimiciv/0.4/).

We utilize OMOP vocabulary to map nonstandard diagnosis codes (primarily ICD from MIMIC) to their standard `SNOMED` equivalent. For medications, we map `NDC` codes to their `RxNorm` equivalent, then to `SNOMED` to ensure that we are capturing the drug ingredient. 

We leverage the hierarchical vocabulary structure to compare indication matchings between:

1\. One-to-one Nonstandard `ICD` to standard `SNOMED` diagnosis mapping,

2\. Horizontal mapping from Standard `SNOMED` to other non-standard `SNOMED`, and

3\. Hierarchical mapping from Standard `SNOMED` to its parents and children vocabulary.

We also leverage the MeSH hierarchy to compare: 

1\. One-to-one MeSH terms 
2\. Tree taversal of MeSH structure

## Data Sources

### Indication Data

- [Drug Central](https://drugcentral.org/): A database that includes formal, FDA-approved indications for each drug, as well as known off-label uses and contraindications. 

- [AACT](https://aact.ctti-clinicaltrials.org/schema): A clinical trial database that contains information on all clinical trials, the intervention being studied (for our purspsoses, the drug) and the conditions that intervention is studied as treating

- [NDF-RT](https://bioportal.bioontology.org/ontologies/NDFRT): Included in the OMOP vocabulary. Developed by the Veterans Health Association, and includes indication relationships within the vocabulary. 

### Brief Description of MIMICIV - Clinical Data Source

MIMICIV is a de-identified, publicly available dataset originated from intensive care units at the Beth Israel Deaconess Medical Center (BIDMC). MIMICIV uses `ICD 9 / 10 (CM)` as their diagnosis codes, and `NDC` codes for medication prescriptions. Each patient is given a unique `subject_id`, and a unique `hadm_id` for each visit. We extract all `ICD` diagnoses and `NDC` medications for each `subject_id` and `hadm_id`, and map them both to `SNOMED` for convenience of vocabulary traversal.

MIMICIV also includes patient demographic and insurance information. We use this information to conduct preliminary implication analysis based on our results.

## Resources Navigation

### Patient_vocab_dict

`off_label/patient_vocab_dict`

Extracted and mapped patient data store in json files. Processing time to extract these is intensive, and so storing them in jsons was most tractable solution for re-accessing original analysis data

### Mapping

`off_label/mapping`

`off_label.py` : Extracts relevant data from MIMICIV publicly available clinical data, which can be accessed at [<https://physionet.org/content/mimiciv/0.4/>](https://physionet.org/content/mimiciv/0.4/){.uri}.\
`vocab_mappings.py` : Python script that defines functions that extract standard to non-standard terms or vice versa, as well as extracting hierarchical concepts. Used to convert from ICD to SNOMED, and to MeSH codes, and the associated hierarchical terms.\

### Results

`off_label/Results`

Finalized summary of Drug Central, AACT, and NDFRT analysis. Includes a powerpoint presentation that summarizes our project. `figures/` include visualizations of our result, and our initial implication analysis.


### Code

`off_label/Code`

`DC_AACT_Results_Aggregation.ipynb`: Loading in jsons of patient data and converting to summary statistics for DrugCentral and AACT matches

`AACT_DC_Ind_Extraction.ipynb`: Pulling all indications listed for drug in question from AACT and DrugCentral databases - interfacing with SQL 

----------

# To Do
* Leverage existing code of extracting all descendants to expand from 1-degree to full hierarchy https://github.com/WengLab-InformaticsResearch/ehr_prevalence/blob/master/ehr_prevalence.py
