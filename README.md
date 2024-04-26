# ehbs_analysis
Code for the paper "Stratifying Risk of Alzheimer’s Disease in Healthy Middle-Aged Individuals with Machine Learning"

The paper aims to study the risk for AD in healthy middle-aged people with abnormal levels of CSF amyloid-beta and Tau, referred as Asymptomatic AD population (figure a). It identifies the CSF peptides which are differentiated between healthy controls and AD subjects (figure b), and uses those peptides to stratify asymptomatic AD (figure c and d) in two datasets (EHBS and ADNI).
<img width="612" alt="image" src="https://github.gatech.edu/storage/user/16378/files/573185fc-0244-4f1d-9109-ef308ab63bf7">

# Peptide selection in EHBS cohort
Code for peptide selection in EHBS and their utility for startifying AsymAD cases is presented in the ipynb file ``peptide_selection_ehbs.ipynb``. The replication of these peptides in the ADNI data is presented in ``ADNI_validation.ipynb``. They both use functions defined in `helper_fn.py`.

# Event-Based Modeling on EHBS and ADNI data
The identified peptides are shown to stratify the AsymAD cases using an alternate modeling approach -the Event Based Model (EBM). The effectiveness of the peptide panel with EBM is shown on both EHBS and ADNI data in the ipynb notebook ``EBM-EHBS-ADNI-dataset.ipynb``. It uses some supporting functions from `analysis_functions.py`

