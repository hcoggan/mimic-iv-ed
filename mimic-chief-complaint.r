library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(testit)
library(httpgd)
hgd()

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(rlang)
library(stringr)

#This file deals with the 'chief complaint data'

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

edstays <- read.csv("edstays_mult_visit_with_vitals.csv")

#convert everything to upper case
edstays$chiefcomplaint <- toupper(edstays$chiefcomplaint) 

#replace some unhelpful acronyms
ha <- which(edstays$chiefcomplaint=="HA")
ah <- which(edstays$chiefcomplaint=="AH")
od <- which(edstays$chiefcomplaint=="OD")
si <- which(edstays$chiefcomplaint=="SI")
st <- which(edstays$chiefcomplaint=="ST")
edstays$chiefcomplaint[ha] <- "H/A"
edstays$chiefcomplaint[ah] <- "A/H"
edstays$chiefcomplaint[od] <- "OVERDOSE"
edstays$chiefcomplaint[si] <- "SUICIDAL"
edstays$chiefcomplaint[st] <- "TACHYCARDIA"



complaints_by_number <- table(edstays$chiefcomplaint)

#there are 55211 different complaints before splitting by comma
#now there are 24090 after splitting by comma (properly complaint-components)
#remove leading spaces, probably resulting from the split
complaints <- edstays$chiefcomplaint
complaints <- unlist(strsplit(complaints,", ")) #comma followed by space
complaints <- unlist(strsplit(complaints,","))
complaints_by_number <- sort(table(complaints), decreasing=TRUE)

write.csv(complaints_by_number, "chief_complaints_by_number.csv")




#only the first 95 complaint-components
#go through and create a dictionary to group together complaints; don't include the key in the associated list


complaint_dict <- list(
    #if a complaint includes this string, we count it as including abdo-pain
    #also include abdominal distention
    "abd_pain" = c("ABDOMINAL PAIN", "ABDO PAIN", "ABD PAIN", "ab PAIN", "EPIGASTRIC PAIN", "ABDOMINAL DISTENTION", "PELVIC PAIN", "RIB PAIN", "RLQ PAIN", "RUQ PAIN"),
    #if a complaint includes chest pain we assume that's enough to count it
    "chest_pain" = c("CHEST PAIN", "CP", "C/P"),
    "transfer" = c("TRANSFER"),
    "fall" = c("FALL"), #this will get "s/p fall" and similar
    "dyspnea" = c("DYSPNEA", "SOB", "SHORTNESS OF BREATH", "SHORT OF BREATH", "WHEEZING", "DIFFICULTY BREATHING", "RESPIRATORY DISTRESS"),
    "headache" = c("HEADACHE", "H/A", "HEAD ACHE", "MIGRAINE"),
    "fever" = c("FEVER"),
    #group together all patients presenting with intoxication/detox
    "intox" = c("ETOH", "SUBSTANCE", "DETOX", "WITHDRAWAL"),
    "overdose" = c("OVERDOSE"),
    #group together nausea, vomiting and diarrhea- N/V will get NVD
    "nvd" = c("N/V", "NAUSEA", "DIARRHEA", "VOMITING"),
    "back_pain" = c("BACK", "LBP"), #this will also get lower back pain, similar
    #group together weakness, dizziness, fainting, fatigue
    "weak_and_dizzy" = c("DIZZINESS", "WEAKNESS", "LIGHTHEADED", "SYNCOPE", "PRESYNCOPE", "FATIGUE", "LETHARGY", "DIZZY", "VERTIGO"),
    "wound" = c("WOUND", "LACERATION"), #this will get WOUND EVAL, WOUND CHECK etc
    #group together people who are presenting for suicidal ideation or a suicide attempt
    "suicidal" = c("SUICIDAL", "SELF HARM", "ATTEMPT"), #almost everyone with 'attempt' in their records has attempted suicide, except one 'failed attempt at self-care'
    "cough" = c("COUGH"),
    #sent to ER following previous exam 
    "abnormal_labs" = c("ABNORMAL", "ABNL"),
    "ili" = c("ILI", "INFLUENZA", "FLU SYMPTOMS", "FLU LIKE", "FLU-LIKE"),
    #confusion/altered mental status
    "confusion" = c("ALTERED MENTAL STATUS", "CONFUSION", "AGITATION"),
    #mental illness, not suicidality or confusion
    "mental_illness" = c("ANXIETY", "DEPRESSION", "PSYCH", "HALLUCINATIONS", "A/H", "DELUSION", "MANIA"),
    #car accidents, inc pedestrians and people hit by bicycles
    "mvc" = c("MVC", "PED STRUCK", "BICYCLE ACCIDENT", "MOTORCYCLE ACCIDENT"),
    "head_injury" = c("HEAD INJURY", "HEAD BLEED", "HEADBLEED", "HEAD LAC"),
    "sore_throat" = c("SORE THROAT"),
    #heart flutters, palpitations
    "palpitations" = c("PALPITATIONS", "TACHYCARDIA", "ABNORMAL EKG", "BRADYCARDIA", "FIBRILLATION"),
    "rectal_pain" = c("RECTAL PAIN"),
    #blood in urine, stool, or vomit, or coughing up blood
    "rectal_bleeding" = c("BRBPR", "RECTAL BLEEDING", "BLEEDING FROM RECTUM", "BLOOD IN STOOL", "BLOODY STOOL", "MELENA", "HEMATURIA", "BLOOD IN URINE", "HEMATEMESIS", "HEMOPTYSIS", "COFFEE GROUND EMESIS"),
    #neck pain/swelling
    "neck_pain" = c("NECK"),
    #seizures
    "seizures" = c("SEIZURE"),
    "rash" = c("RASH"),
    #painful or difficult urination
    "dysuria" = c("DYSURIA", "URINARY RETENTION", "URINARY FREQUENCY"),
    #this will get both left and right flank pain
    "flank_pain" = c("FLANK"),
    "walking_problems" = c("UNSTEADY GAIT", "UNABLE TO AMBULATE", "DIFFICULTY AMBULATING", "INABILITY TO AMBULATE", "DIFF AMBULATING", "CAN'T WALK", "CANT WALK", "UNABLE TO WALK", "VAGINAL SPOTTING"),
    "hypertension" = c("HYPERTENSION", "HIGH BLOOD PRESSURE"),
    #vaginal bleeding/pain
    "vaginal_bleeding" = c("VAGINAL BLEEDING", "VAG BLEEDING", "VAG BLEED", "VAGINAL BLEED", "VAGINAL PAIN", "VAGINAL SPOTTING",  "VAG SPOTTING"),
    #group together blood sufar issues
    "blood_sugar" = c("HYPERGLYCEMIA", "HYPOGLYCEMIA", "BLOOD SUGAR"),
    "allergy" = c("ALLERGIC", "ALLERGY"),
    "pregnant" = c("PREG", "LABOR"), #possibly this will get people 'requesting pregnancy test', but there are very very few of those
    "knee_pain" = c("KNEE"), #again will get people with pain on both sides
    "assault" = c("ASSAULT"),
    "low_blood_pressure" = c("HYPOTENSION", "LOW BLOOD PRESSURE"),
    #group together ankle, toe and foot pain
    "foot_pain" = c("FOOT", "ANKLE", "TOE"),
    #any vision issues
    "vision_issues" = c("VISION", "SEEING FLOATERS", "BLIND", "VISUAL"),
    #in general, group pain and swelling
    "arm_pain" = c("ARM", "ELBOW"),
    "shoulder_pain" = c("SHOULDER"),
    "hip_pain" = c("HIP"),
    "abscess" = c("ABCESS"),
    "asthma" = c("ASTHMA"),
    "nosebleed" = c("EPISTAXIS"),
    "hypoxia" = c("HYPOXIA"),
    #brain hemorrhages
    "hematoma" = c("ICH", "SDH", "SAH"),
    #group together pain and swelling
    "leg_pain_or_swelling" =c("LEG", "CALF"),
    "eye_pain" = c("EYE PAIN"), #multiple categories
    "ear_pain" = c("EAR PAIN"),
    "numbness" = c("NUMBNESS"),
    "hand_finger_wrist" = c("HAND", "FINGER", "WRIST")
)


print(head(grepl("ABD PAIN", edstays$chiefcomplaint)))

#also create a specialised file just for complaints
patient_complaints <- data.frame(stay_id=edstays$stay_id)

#now go through and assign them
categs_of_complaint <- names(complaint_dict)
print(length(categs_of_complaint))

for (name in categs_of_complaint) {
    edstays[[name]] <- c(rep(FALSE, nrow(edstays)))
    strings <- complaint_dict[[name]]

    for (string in strings) {
        edstays[[name]] <- edstays[[name]] | grepl(string, edstays$chiefcomplaint)
    }

    patient_complaints[[name]] <- edstays[[name]]

}

#drop chief complaint column
edstays <- subset(edstays, select = -c(chiefcomplaint))

#how many patients have at least one recorded complaint?
fraction_described <-  mean(rowSums(edstays[categs_of_complaint]) > 0)


print("Fraction of patients with legible complaint:")
print(fraction_described)

print(head(edstays))

write.csv(edstays, "edstays_with_binarised_complaints.csv")
write.csv(patient_complaints, "binarised_complaints_alone.csv")


