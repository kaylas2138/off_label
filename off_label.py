import numpy as np
import pandas as pd
import os
# get your current working directory
os.getcwd()
# set your working directory to the directory in which you have all your mimic data files
os.chdir('/Users/yc3972/Desktop/DBMI/Courses/GD1_Fall/G4003 Symbolic Methods/Project/off_label')

# Reading in datasets
admissions = pd.read_csv("./clinical_data/physionet.org/files/mimiciv/1.0/core/admissions.csv", header = 0, sep = ',', engine = 'c')
admissions['admittime'] = pd.to_datetime(admissions['admittime'])
admissions['dischtime'] = pd.to_datetime(admissions['dischtime'])
admissions['insurance'] = admissions['insurance'].astype('|S') 
admissions['duration'] = admissions['dischtime'] - admissions['admittime']
prescriptions = pd.read_csv("./clinical_data/physionet.org/files/mimiciv/1.0/hosp/prescriptions.csv",header =0,sep = ',', 
                            dtype={'drug':'str','ndc': 'str'})
diagnosis = pd.read_csv("./clinical_data/physionet.org/files/mimiciv/1.0/hosp/diagnoses_icd.csv", sep = ',', header = 0)
diagnosis_mapping = pd.read_csv("./clinical_data/physionet.org/files/mimiciv/1.0/hosp/d_icd_diagnoses.csv", sep = ',', header =0)
diagnosis2 = pd.merge(diagnosis, diagnosis_mapping, on = ['icd_code', 'icd_version'])

# patient IDs and visit IDs of patients who have been prescribed with given drug
def patient_list(drug):
    patients = prescriptions.loc[prescriptions.drug.str.contains(drug, case = False, na = False, regex = False)][['subject_id', 'hadm_id']].reset_index()
    return patients

# list of all hospital visits from a specific patient (that may not have included drug of query), for generalizability
def unique_visits(patient_id):
    visits = admissions.loc[admissions.subject_id == patient_id]['hadm_id'].unique()
    return visits

# list of all the diagnosis given for that visit
def diagnosis_list(patient_id, visit_id):
    all_diagnosis = diagnosis2.loc[diagnosis2.subject_id == patient_id]
    visit_diagnosis = all_diagnosis.loc[all_diagnosis.hadm_id == visit_id].sort_values(by = 'seq_num').reset_index()
    return visit_diagnosis

# list of all the drugs given for that visit
def drugs_list(patient_id, visit_id):
    all_drugs = prescriptions.loc[prescriptions.subject_id == patient_id]
    visit_drugs = all_drugs.loc[all_drugs.hadm_id == visit_id].reset_index()
    return visit_drugs


# For the sake of simplicity, we will be investigating one of the patients who have been given this drug
def returning_lists(drug, i):
    patients_given = patient_list(drug)
    if i > patients_given.shape[0]:
        raise ValueError("Given index is bigger than size of patients list of the drug")
    query_patient = patients_given.subject_id[i]
    query_visit = patients_given.hadm_id[i]
    print("Searching for Patient "+ str(query_patient) + " who has been given " + str(drug) + " , for visit number " + str(query_visit) + ".")
    diagnoses = diagnosis_list(query_patient, query_visit)
    drugs = drugs_list(query_patient, query_visit)
    return diagnoses, drugs

# Understand descriptions of the patient: gender and insurance information 
def patient_details(subject_id_list):
    details = admissions[admissions.subject_id.isin(subject_id_list)].reset_index()
    return details

