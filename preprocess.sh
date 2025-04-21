#!/bin/bash
set -e

# 1. "edstays-extract-time-of-day-and-previous-visits.r": extracts time spent in ER, time of arrival, year group associated with patient, and visit history.
# 2. "edstays_triage_data_cleaning_mult_visits.r": cleans and categorises vital signs, pain scores, and demographic variables.
# 3. "mimic-chief-complaint.r": translates the free-text "chief complaint" field to a series of 54 binary indicators associated with each visit.
# 4. "link-with-mimic-iv.r": obtains patient age at each visit from the MIMIC-IV dataset, and translates ICD codes associated with each stay to a Charlson score. Also obtains insurance, marital, and language info for visits resulting in an admission.
# 5. "edstays-triage-hists.r": produces histograms showing the distribution of variables across the cohort, allowing at-a-glance identification of reference groups.
# 6. "mimic_mult_visit_with_vitals_binary_recode.r": reformat data on all visits to ER to a series of binary indicators.
# 7. "adm_edstays_binary_recode.r": reformat data on all visits to ER *resulting in an admission* to a series of binary indicators, including additional covariates.
# 8. "edstays_add_vitalsigns_during_stay.r": attach data on vital signs from readings taken during a patient's stay, cleaned, categorised, and reformat as binary indicators.
mkdir -p plots
steps=(
#  "edstays-extract-time-of-day-and-previous-visits.r"
#  "edstays_triage_data_cleaning_mult_visits.r"
#  "mimic-chief-complaint.r"
#  "link-with-mimic-iv.r"
#  "edstays_triage_hists.r"
#  "mimic_mult_visit_with_vitals_binary_recode.r"
#  "adm_edstays_binary_recode.r"
 "edstays_add_vitalsigns_during_stay.r"
 )


for step in ${steps[@]} ; do 
    echo "########################################"
    echo ${step} ...
    echo "####################"
    Rscript ${step}
    echo "########################################"
done