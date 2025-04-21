library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(testit)
library(httpgd)
library(broom)
library(stringr)
library(patchwork)
library(forcats)
library(pROC)
library(comorbidity)
hgd()

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(janitor)

#This file retrieves age, language, insurance and marital status from MIMIC database
#we also link ICD code from MIMIC

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

edstays <- read.csv("edstays_with_binarised_complaints.csv")
#from MIMIC data
patients <- read.csv("mimic-iv-3.1/hosp/patients.csv.gz")
admissions <- read.csv("mimic-iv-3.1/hosp/admissions.csv.gz")
diagnosis <- read.csv("mimic-iv-ed-2.2/ed/diagnosis.csv")


#triple check this: how many patients included in admissions are included in edstays?
adm_subj_ids <- unique(admissions$subject_id)
eds_subj_ids <- unique(edstays$subject_id)

print("Number of patients included in edstays")
print(length(eds_subj_ids))


print("Fraction of patients included in admissions table")
with_insurance_info <- intersect(adm_subj_ids, eds_subj_ids)
print(length(with_insurance_info)/length(eds_subj_ids))


stays_with_ins_info <- which(edstays$subject_id %in% with_insurance_info)
adm_decision <- edstays$admitted[stays_with_ins_info]

print("Of stays with a subject_id linked in the admissions table,")
print(sum(adm_decision)/length(adm_decision))
print("result in an admission")

print("relative to ")
print(sum(edstays$admitted)/nrow(edstays))
print("of all stays")


#recode everything written as TRUE/FALSE 
edstays[edstays==TRUE] <- 1
edstays[edstays==FALSE] <- 0
 
#reduce patient data down to the data we want (just anchor-age, which should be considered only in bins larger than 3
#due to the size of the anchor-bins, and caps out at 91)

#probably 10 year bins
#anyway, retrieve age and year
patients <- patients %>% select(subject_id, anchor_age, anchor_year, anchor_year_group, dod)


#only one age per patient
assert(length(unique(patients$subject_id))==nrow(patients))
print("Assertion passed")

edstays <- edstays %>% inner_join(patients, by="subject_id")
#assume everyone is born on 1 Jan
under_91 <- which(edstays$anchor_age <= 90)
over_90 <- which(edstays$anchor_age == 91) #everyone over 89 in anchor year has age set to 91

edstays$age <- 91 #maintaining anchor-age settings by default

edstays$age[under_91] <- edstays$anchor_age[under_91] + edstays$year_of_arrival[under_91] - edstays$anchor_year[under_91]

assert(min(edstays$anchor_age)>=18)
assert(min(edstays$age)>=18)
print("No children")

print(head(edstays))
#how many patients have an age recorded?
age_na <- sum(is.na(edstays$age))/nrow(edstays)
print("fraction with no recorded age")
print(age_na)

#calculate precise age from admitted year and anchor year

#group each age set into rough decades
age_labels <- c("18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+")
lower_limits <- c(18, 30, 40, 50, 60, 70, 80, 90)
upper_limits <- c(29, 39, 49, 59, 69, 79, 89, 150)

#create bins
edstays$age_group <- c(rep(0, nrow(edstays)))


for (i in 1:length(age_labels)) {
    edstays$age_group[which(between(edstays$age, lower_limits[i], upper_limits[i]))] <- age_labels[i]
}

edstays$intime_as_date <- as.POSIXct(edstays$intime_as_date, format = "%Y-%m-%d")


#now link Charlson score
#first combine version with code- take MOST RELEVANT CODE only
diagnosis <- diagnosis %>% filter(seq_num==1)
diagnosis$charlson_score <- NA
#differentiate between ICD-9 and ICD-10
v9 <- which(diagnosis$icd_version == 9)
v10 <- which(diagnosis$icd_version == 10)

#assign comorbidities with hierarchy
com9 <- comorbidity(diagnosis[v9,], id='stay_id', code='icd_code', map="charlson_icd9_quan", assign0=TRUE)
com10 <- comorbidity(diagnosis[v10,], id='stay_id', code='icd_code', map="charlson_icd10_quan", assign0=TRUE)

diagnosis$charlson_score[v9] <- comorbidity::score(com9, w="quan", assign0=TRUE)
diagnosis$charlson_score[v10] <- comorbidity::score(com10, w="quan", assign0=TRUE)

diagnosis <- diagnosis %>% select(stay_id, charlson_score)

print("Visits before ICD join")
print(nrow(edstays))
edstays <- edstays %>% inner_join(diagnosis, by='stay_id')
print("Visits after ICD join")
print(nrow(edstays))

#link charlson values
print("Charlson table")
print(table(edstays$charlson_score))


#now link for admitted patients
adm_edstays <- edstays %>% filter(admitted==1)

admissions <- admissions %>% select(hadm_id, admittime, admission_location, insurance, language, marital_status)

print("before linkage")
print(nrow(adm_edstays)) #138564 visits
adm_edstays <- adm_edstays %>% inner_join(admissions, by='hadm_id')
print("after linkage")
print(nrow(adm_edstays)) #138219 visits



print(head(adm_edstays))

#check all were admitted from ED
print(table(adm_edstays$admission_location))
#ok, they weren't; check with Bill

#fill N/As instead of blank spaces
blank <- which(adm_edstays$marital_status == "")
adm_edstays$marital_status[blank] <- 'N/A'

blank <- which(adm_edstays$language== "")
adm_edstays$language[blank] <- 'N/A'

blank <- which(adm_edstays$insurance == "")
adm_edstays$insurance[blank] <- 'N/A'

#to condense language- record whether english is primary language
adm_edstays$primary_english <- rep(0, nrow(adm_edstays))
adm_edstays$primary_english[which(adm_edstays$language=="English")] <- 1

#condense insurance status into four groups- medicare, medicaid, private, other
adm_edstays$insurance_group <- adm_edstays$insurance
to_condense <- which((adm_edstays$insurance_group=="Other") | (adm_edstays$insurance_group=="N/A") | (adm_edstays$insurance_group=="No charge"))
adm_edstays$insurance_group[to_condense] <- "Other"





#drop unnecessary columns
edstays <- subset(edstays, select = -c(anchor_year, anchor_age))

#copy this for linking on date of death
edstays_for_dod <- edstays
adm_edstays_for_dod <- adm_edstays


edstays <- subset(edstays, select = -c(dod, intime_as_date))
adm_edstays <- subset(adm_edstays, select = -c(admitted, discharged, admittime, anchor_year, anchor_age, dod, intime_as_date))

write.csv(edstays, "edstays_with_binarised_complaints_and_age.csv")
write.csv(adm_edstays, "adm_edstays.csv")

#Nothing below this matters for the sake of AMIA draft


#OK! Now filter out patients with any visit history and record time of death if available
edstays_for_dod <- edstays_for_dod %>% 
    filter(discharged==1) %>%
    group_by(subject_id) %>%
    arrange(intime) %>%
    slice_head(n=1) #keep only the first visit



edstays_for_dod$dod <- as.POSIXct(edstays_for_dod$dod, format = "%Y-%m-%d")
edstays_for_dod$time_to_death <- difftime(edstays_for_dod$dod, edstays_for_dod$intime_as_date, units="days")

assert(length(which(!is.na(edstays_for_dod$time_to_death)))>0)

edstays_for_dod <- subset(edstays_for_dod, select = -c(dod, intime_as_date))

write.csv(edstays_for_dod, "edstays_first_visit_discharged_with_linked_death.csv")


#OK! Now for ADMITTED patients filter out patients with any visit history and record time of death if available
adm_edstays_for_dod <- adm_edstays_for_dod %>% 
    group_by(subject_id) %>%
    arrange(intime) %>%
    slice_head(n=1) #keep only the first visit


adm_edstays_for_dod$dod <- as.POSIXct(adm_edstays_for_dod$dod, format = "%Y-%m-%d")
adm_edstays_for_dod$time_to_death <- difftime(adm_edstays_for_dod$dod, adm_edstays_for_dod$intime_as_date, units="days")

adm_edstays_for_dod <- subset(adm_edstays_for_dod, select = -c(admitted, discharged, admittime, anchor_year, anchor_age, dod, intime_as_date))

write.csv(adm_edstays_for_dod, "adm_edstays_first_visit_with_linked_death.csv")
