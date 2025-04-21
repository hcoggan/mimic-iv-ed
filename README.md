# mimic-iv-ed

A repository of code used to analyse the MIMIC-IV-ED dataset

If you find this code useful, please cite:

Helena Coggan, Pradip Chaudhari, Yuval Barak-Corren, Andrew M. Fine, Ben Y. Reis, Jaya Aysola, William G. La Cava.
Demographic Factors Associated with Triage Acuity, Admission and Length of Stay During Adult Emergency Department Visits.
https://arxiv.org/abs/2503.22781

# Study design

We conducted a single-center retrospective cross-sectional study of ED visits by adults to BIDMC, recorded in the MIMIC-IV-ED dataset (in table `edstays`). Visits were excluded if they did not have a valid recorded `in-time` (time of arrival to the ED) or `out-time` (time at which they left the ED); to be valid, an in-time had to precede the recorded out-time. We obtained vital signs taken at triage (from table `triage`) and excluded visits for which such data was not available, or where vitals fell outside a reasonable range (e.g. a temperature of 200 degrees Farenheit). We obtained vital signs recorded during a patient's stay (from table `vitalsigns`), but included stays without such information, as many patients may be admitted too quickly to allow repeated vital assessments. 

Age was obtained to the nearest year by linkage to MIMIC-IV (from table `patients`), assuming a birthday of January 1st for each patient. The group of years associated with each patient's care was obtained from the same table, and could precede 2011; as only one year group could be associated with each patient regardless of the time frame of their visits, we took this as a rough ordering of the time periods in which different patients engaged with the hospital system, rather than a precise temporal identifier. We took ICD codes relating to the diagnosis arising from each ED stay from table `diagnosis`, considering only the most relevant code, and converted these to Charlson comorbidity scores (a measure of severity) using the R package `comorbidity` (Gasparini 2018). Visits were excluded if no ICD code was available; linked ages were observed for all patients.

Patient-reported pain, a free-text field, was converted to a score between 0 and 10, or recorded as "unable to answer" or "non-numeric response" if such a conversion was impossible. Chief complaints were converted to 54 manually curated binary indicators covering the most frequent reasons given for a visit to the ER, such that 88% of visits were described by at least one variable.

Visits resulting in an admission to hospital were assigned a specific linkage ID which allowed us to obtain data on the patient's insurance status, marital status, and primary language from the MIMIC-IV dataset. These data were not available for visits which did not result in an admission. As a result, we were only able to control for these variables when examining demographic disparities in LOS amongst those admitted to the hospital.

# Usage 


## Requirements

- R version 4.4.2 or higher
- `packages.txt` contains the R libraries that are needed. To install them, run `Rscript install_packages.r`. 
- [MIMIC-IV-ED dataset (v 2.2) from Physionet](https://physionet.org/content/mimic-iv-ed/2.2/), mapped to the directory `mimic-iv-ed` in the repo root. 
- [MIMIC-IV dataset (v 3.1) from Physionet](https://physionet.org/content/mimiciv/3.1/), mapped to the directory `mimic-iv-3.1` in the repo root. 

## Preprocessing

These files, run in order, reproduce the analysis we've used to analyse the MIMIC-IV-ED dataset.
To run them all, run `preprocess.sh` (from a *nix system). 

1. `edstays-extract-time-of-day-and-previous-visits.r`: extracts time spent in ER, time of arrival, year group associated with patient, and visit history.
2. `edstays_triage_data_cleaning_mult_visits.r`: cleans and categorises vital signs, pain scores, and demographic variables.
3. `mimic-chief-complaint.r`: translates the free-text `chief complaint` field to a series of 54 binary indicators associated with each visit.
4. `link-with-mimic-iv.r`: obtains patient age at each visit from the MIMIC-IV dataset, and translates ICD codes associated with each stay to a Charlson score. Also obtains insurance, marital, and language info for visits resulting in an admission.
5. `edstays-triage-hists.r`: produces histograms showing the distribution of variables across the cohort, allowing at-a-glance identification of reference groups.
6. `mimic_mult_visit_with_vitals_binary_recode.r`: reformat data on all visits to ER to a series of binary indicators.
7. `adm_edstays_binary_recode.r`: reformat data on all visits to ER *resulting in an admission* to a series of binary indicators, including additional covariates.
8. `edstays_add_vitalsigns_during_stay.r`: attach data on vital signs from readings taken during a patient's stay, cleaned, categorised, and reformat as binary indicators.

## Propensity score and Odds Ratio Calculations
All other files (starting "ps2") are named according to their exposure and outcome. 
These fit propensity scores for a given exposure variable, calculate odds ratios using two PS methodologies (covariate adjustment and matching with McNemar's test), save them, and produce OR plots (for categorical exposures).
To run all of these files, run `analysis.sh`. 

# Acknowledgments

This work is supported in part by the National Library Of Medicine of the National Institutes of Health under Award Number R01LM014300. 
The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.