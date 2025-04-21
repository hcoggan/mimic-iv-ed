#!/bin/bash
# set -e

# run all the ps2 r files
# steps=$(ls ps2_*.r)
steps=(
"ps2_age_acuity.r"
"ps2_age_admission.r"
"ps2_age_adm_waittime.r"
"ps2_age_discharged_waittime.r"
"ps2_gender_acuity.r"
"ps2_gender_admission.r"
"ps2_gender_waittime_adm.r"
"ps2_gender_waittime_discharged.r"
"ps2_race_acuity.r"
"ps2_race_admissions.r"
"ps2_race_age_acuity.r"
"ps2_race_age_admission.r"
"ps2_race_age_adm_waittime.r" #matchit error
"ps2_race_age_discharged_waittimes.r"
"ps2_race_waittime_adm.r"
"ps2_race_waittime_discharged.r"
"ps2_insurance_waittime_adm.r"
"ps2_intersection_adm_waittime.r"
"ps2_intersection_discharged_waittime.r"
"ps2_language_waittime_adm.r"
"ps2_gender_by_race_acuity.r"
"ps2_gender_by_race_adm_discharged.r"
"ps2_gender_by_race_admission.r"
"ps2_gender_by_race_adm_waittimes.r"
"ps2_gender_by_race_waittimes.r" #matchit error
"ps2_rg_interaction_acuity.r"
"ps2_rg_intersection_acuity.r"
"ps2_rg_intersection_admission.r"
)
echo ${steps}

N=8

for step in ${steps[@]} ; do 
        echo "########################################"
        echo  ${step} ...
        echo "####################"
        Rscript ${step}
        echo "########################################"
done
