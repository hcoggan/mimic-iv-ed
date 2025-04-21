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

#This file combines triage and edstays files into a single CSV, cleans up 'vitals'
#and extracts time spent in ER
#takes patients from MULTIPLE visits
#uses same filtration as La Cava github

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")


#assume that we have run edstays-extract-time-of-day-and-prev visits
#so we have already extracted time of arrival, and visits in last several weeks/months
edstays <- read.csv("edstays_mult_visit_with_timeframes_without_vitals.csv")
triage <-  read.csv("mimic-iv-ed-2.2/ed/triage.csv")

print("dataframes loaded")


#first, time spent in ED, measured in minutes
time_in_ed <- edstays$outtime - edstays$intime
edstays$time_in_ed <- time_in_ed

print('calculated time in ED')

zero_vector <- c(rep(0, nrow(edstays)))

#create a new column just recording if admission happened
#also recording if patient was discharged home
edstays$admitted <- zero_vector
edstays$discharged <- zero_vector

edstays$admitted[which(edstays$disposition=='ADMITTED')] <- 1
edstays$discharged[which(edstays$disposition=='HOME')] <- 1


print('calculated discharge/admission')

#output unique racial coding
race_coding <- unique(edstays$race)

#recode races following the lead of the La Cava github:
 race_recode <- c("WHITE" = "White", "ASIAN" = "Asian", "OTHER" = "Other", "ASIAN - SOUTH EAST ASIAN" = "Asian",
     "UNKNOWN" = "Unknown", "BLACK/AFRICAN AMERICAN" = "Black/African-American", "ASIAN - CHINESE" = "Asian", "BLACK/CAPE VERDEAN" = "Black/African-American", 
     "PORTUGUESE" = "Other", "MULTIPLE RACE/ETHNICITY" = "Unknown", "UNABLE TO OBTAIN" = "Unknown", "WHITE - OTHER EUROPEAN" = "White", 
     "WHITE - RUSSIAN" = "White", "HISPANIC/LATINO - DOMINICAN" = "Hispanic/Latino", "HISPANIC/LATINO - SALVADORAN" = "Hispanic/Latino", "HISPANIC/LATINO - PUERTO RICAN" = "Hispanic/Latino",
      "WHITE - BRAZILIAN" = "White", "HISPANIC/LATINO - GUATEMALAN" = "Hispanic/Latino", "BLACK/CARIBBEAN ISLAND" = "Black/African-American", "HISPANIC OR LATINO" = "Hispanic/Latino", "BLACK/AFRICAN" = "Black/African-American", 
      "ASIAN - KOREAN" = "Asian", "WHITE - EASTERN EUROPEAN" = "White", "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" = "Native Hawaiian/Other Pacific Islander", "HISPANIC/LATINO - MEXICAN" = "Hispanic/Latino", 
      "PATIENT DECLINED TO ANSWER" = "Unknown", "HISPANIC/LATINO - CUBAN" = "Hispanic/Latino", "HISPANIC/LATINO - COLUMBIAN" = "Hispanic/Latino", "ASIAN - ASIAN INDIAN" = "Asian", "AMERICAN INDIAN/ALASKA NATIVE" = "American Indian/Alaska Native", 
      "SOUTH AMERICAN" = "Other", "HISPANIC/LATINO - HONDURAN" = "Hispanic/Latino", "HISPANIC/LATINO - CENTRAL AMERICAN" = "Hispanic/Latino")


#recode race 
edstays$race <- dplyr::recode(edstays$race, !!! race_recode)

#output unique arrival coding
arrival_coding <- unique(edstays$arrival_transport)
for (arrival in arrival_coding) {
    print(arrival)
    print(length(which(edstays$arrival_transport == arrival))) #how many people were recorded as arriving by each method?
}

print("we have recoded race and arrival")
print(colnames(triage))


#drop subject_id from triage (we have in it edstays)
triage <- triage %>% select(stay_id, temperature, heartrate, resprate, o2sat, sbp, dbp, pain, acuity, chiefcomplaint)

edstays <- edstays %>% inner_join(triage, by="stay_id")

assert(nrow(edstays)==length(unique(edstays$stay_id)))
print("Duplicate rows not reported per stay")


print("unique pain values")
print(length(unique(edstays$pain)))


#first, if something is clearly intended as a number between 0 and 10, keep it
valid_pain_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")


#finally, condense every integer up to 20 down to 10
indices <- which((!is.na(edstays$pain)) & (as.numeric(edstays$pain) < 21) & (as.numeric(edstays$pain) > 10))
edstays$pain <- replace(edstays$pain, indices, "10")

for (pain in valid_pain_levels) {
    #find those with random punctuation either side
    intended_level <- which((edstays$pain == paste(" ", pain, sep="")) |
         (edstays$pain == paste(pain, " ", sep="")) | (edstays$pain == paste("  ", pain, sep="")) | 
         (edstays$pain == paste(pain, "  ", sep="")) | (edstays$pain == paste(pain, ".", sep="")) | 
         (edstays$pain == paste(".", pain, sep="")) | (edstays$pain == paste(pain, "-", sep="")) | 
         (edstays$pain == paste("-", pain, sep="")) | (edstays$pain == paste(pain, "/10", sep="")) |
         (edstays$pain == paste(pain, ",", sep="")) | (edstays$pain == paste("0", pain, sep="")) | (edstays$pain == paste(",", pain, sep="")))
    edstays$pain <- replace(edstays$pain, intended_level, pain)
}

#first condense intermediate or out of range values
for (pain in valid_pain_levels) {
    for (pain_ in valid_pain_levels) {
        p <- as.numeric(pain)
        p_ <- as.numeric(pain_)
        #find ranges (i.e. '6-7')
        rep_value <- (p+p_)/2
        low_int <- floor(rep_value) #lowest integer
        indices <- which(edstays$pain == paste(pain, pain_, sep="-"))
        if (!(rep_value == low_int)) { #if the middle of this range is not an integer
            prob_higher_integer <- rep_value - low_int #i.e. 3.8 is assigned to 4 80% of the time, and to 3 20% of the time
            for (index in indices) {
                if (runif(1) < prob_higher_integer) { #assign the lower or higher integer probabilistically 
                    edstays$pain[index] <- as.character(low_int + 1)
                } else {
                    edstays$pain[index] <- as.character(low_int)
                }
            }
        } else { #middle of this range is an integer
                for (index in indices) {
                    edstays$pain[index] <- as.character(low_int)
                }
        }
        #now find values between 20 and 100
        val <- 10*p+p_
        if (val > 20) {
            val_char <- as.character(val)
            indices <- which(edstays$pain == val_char)
            #replace it assuming it was meant to be between 1 and 100
            rep_value <- val/10 
            low_int <- floor(rep_value) #lowest integer
            if (!(rep_value == low_int)) { #if the middle of this range is not an integer
                prob_higher_integer <- rep_value - low_int #i.e. 3.8 is assigned to 4 80% of the time, and to 3 20% of the time
                for (index in indices) {
                    if (runif(1) < prob_higher_integer) { #assign the lower or higher integer probabilistically 
                        edstays$pain[index] <- as.character(low_int + 1)
                    } else {
                        edstays$pain[index] <- as.character(low_int)
                    }
                }
            } else { #middle of this range is an integer
                for (index in indices) {
                    edstays$pain[index] <- as.character(low_int)
                }
        }
        }
        #now find intermediate values e.g. 6.5
        val_char <- paste(pain, pain_, sep=".")
        rep_value <- as.numeric(val_char)
        low_int <- floor(rep_value) #integer below it
        indices <- which(edstays$pain == val_char)
        prob_higher_integer <- rep_value - low_int #i.e. 3.8 is assigned to 4 80% of the time, and to 3 20% of the time
        for (index in indices) {
            if (runif(1) < prob_higher_integer) { #assign the lower or higher integer probabilistically 
                edstays$pain[index] <- as.character(low_int + 1)
            } else {
                edstays$pain[index] <- as.character(low_int)
            }
        }

    }   
}



#keep anything that clearly means "unable to answer" (including 'critical')

unable_synonyms <- c("u", "U", "U.A", "u/'a", "u/a", "u/a ", "u/a.  ",
            "ua", "uA", "Ua", "UA", "uable", "uanable", "uanble", "uanble to answer",
            "uanlbe", "uable", "uanable", "uanble", "uanble to answer", 
            "uanlbe", "unabale", "unabe", "unabke", "unabkle", "unabl e", 
            "unabl3", "unable", "Unable", "UNABLE", "unable ", "Unable ",
            "UUTA", "utta", "uts", "ut", "urta",
            "unable 97/", "unable critical", "unable to assess", "unable to obtain",
            "unable to scale", "unable to score", "unable trauma", "unable.",
            "unablwe", "unalbe", "unale", "unale ", "unbale", "unbel", "unble", "uta",
            "uTA", "UTa", "UTA", "uta ", "UTA ", "critical", "CRITICAL", "Critical", "c",
            "crit", "criot", "crit ", "critcal", "critical ", "critial ", "critica;",
            "critical ", "critical; ", "criticial", "criut", "crtiical", "crtitcal")

unable_indices <- which(edstays$pain %in% unable_synonyms)
edstays$pain <- replace(edstays$pain, unable_indices, "UTA")
                       

#print(table(edstays$pain))

print("unique pain values after filtration")
print(length(unique(edstays$pain)))

#now there are 620

#final values
valid_pain_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "UTA")

#condense everything else down to 'Other'
other_indices <- which(!(edstays$pain %in% valid_pain_levels))
edstays$pain <- replace(edstays$pain, other_indices, "non-numeric")

#now condense numeric response into none, mild, moderate, severe
none <- which(edstays$pain=="0")
mild <- which((edstays$pain=="1") | (edstays$pain=="2") | (edstays$pain=="3"))
moderate <- which((edstays$pain=="4") | (edstays$pain=="5") | (edstays$pain=="6") | (edstays$pain=="7"))
severe <- which((edstays$pain=="8") | (edstays$pain=="9") | (edstays$pain=="10"))

edstays$pain <- replace(edstays$pain, none, "none")
edstays$pain <- replace(edstays$pain, mild, "mild")
edstays$pain <- replace(edstays$pain, moderate, "moderate")
edstays$pain <- replace(edstays$pain, severe, "severe")

#condense acuity 
very_high <- which(edstays$acuity=="1")
high <- which(edstays$acuity=="2")
moderate <- which((edstays$acuity=="3")) 
low <- which((edstays$acuity=="4") | (edstays$acuity=="5"))
unassigned <- which(is.na(edstays$acuity))

edstays$acuity <- replace(edstays$acuity, very_high, "very_high")
edstays$acuity <- replace(edstays$acuity, high, "high")
edstays$acuity <- replace(edstays$acuity, moderate, "moderate")
edstays$acuity <- replace(edstays$acuity, low, "low")
edstays$acuity <- replace(edstays$acuity, unassigned, "unassigned")

#filter temperature, resp rate, blood pressure on La Cava Github


edstays <- edstays %>% filter(between(temperature, 95, 105) & between(resprate, 2, 200)
        & between(o2sat, 50, 100) & between(sbp, 30, 400) & between(dbp,30, 300)
     & between(heartrate, 30, 300))


#SAVE DATA BEFORE BINARISATIONS

write.csv(edstays, "edstays_cleaned_continuous.csv")

#for AMIA paper: print out dispositions

print(table(edstays$disposition))

#define fever per Andrew Fine and Pradip Chaudhari
fever_grade <- c(rep("none", nrow(edstays)))
fever <- which(edstays$temperature > 100.4)
not_known <- which(is.na(edstays$temperature))
fever_grade[fever] <- "fever"
fever_grade[not_known] <- "not_known"

edstays$fever_grade <- fever_grade

#define O2 level per Wold 2015 (and Jaya Aysola)
o2_level <- c(rep("normal", nrow(edstays)))
low <- which(edstays$o2sat <= 95.0 & edstays$o2sat > 85)
very_low <- which(edstays$o2sat <= 85)
not_known <- which(is.na(edstays$o2sat))
#very_low <- which(edstays$o2sat <= 92)


o2_level <- replace(o2_level, low, "85_to_95")
o2_level <- replace(o2_level, very_low, "below_85")
o2_level <- replace(o2_level, not_known, "not_known")
#o2_level <- replace(o2_level, very_low, "very low")

edstays$o2_level <- o2_level


#define resp level per NIH
resp_level <- c(rep("normal", nrow(edstays)))
low <- which(edstays$resprate < 12)
high <- which(edstays$resprate > 20)
not_known <- which(is.na(edstays$resp_level))


resp_level <- replace(resp_level, low, "low")
resp_level <- replace(resp_level, high, "high")
resp_level <- replace(resp_level, not_known, "not_known")


edstays$resp_level <- resp_level

#also heart rate levels, per Nursing Skills textbook
heart_level <- c(rep("normal", nrow(edstays)))
low <- which(edstays$heartrate < 60)
high <- which(edstays$heartrate > 100 & edstays$heartrate <= 150)
very_high <- which(edstays$heartrate > 150)
not_known <- which(is.na(edstays$heartrate))


heart_level <- replace(heart_level, low, "low")
heart_level <- replace(heart_level, high, "high")
heart_level <- replace(heart_level, very_high, "very_high")
heart_level <- replace(heart_level, not_known, "not_known")


edstays$heart_level <- heart_level

#also blood pressure levels, per Johns Hopkins

bp_level <- c(rep("normal", nrow(edstays)))
low <- which((edstays$sbp < 90) & (edstays$dbp < 60))
elevated <- which((edstays$sbp >= 120) & (edstays$sbp < 130) & (edstays$dbp < 80))
stage_1_high <- which(((edstays$sbp >= 130) & (edstays$sbp < 140)) | ((edstays$dbp >= 80) & (edstays$dbp < 90)))
stage_2_high <- which(((edstays$sbp >= 140) & (edstays$sbp < 180)) | (edstays$dbp >= 90)) #alter definition to make sure unstable is the right category
unstable <- which(edstays$sbp >= 180)
not_known <- which(is.na(edstays$sbp) | is.na(edstays$dbp))

bp_level <- replace(bp_level, low, "low")
bp_level <- replace(bp_level, elevated, "elevated")
bp_level <- replace(bp_level, stage_1_high, "stage_1_high")
bp_level <- replace(bp_level, stage_2_high, "stage_2_high") #140/90 is stage 2 hypertension threshold-- if something 'could also be' stage 1 by above definition, we assign it stage 2
bp_level <- replace(bp_level, unstable, "unstable") 
bp_level <- replace(bp_level, not_known, "not_known")



#check that everyone not caught by these does in fact have normal blood pressure 
classified_normal <- which(bp_level == "normal")
actually_normal <- which((edstays$sbp < 120) & (edstays$dbp < 80))
print("Set difference") #we want this to be 0
print(length(setdiff(classified_normal, actually_normal)))
print(length)

edstays$bp_level <- bp_level


print(head(edstays))
print("number of visits")
print(nrow(edstays))

#do not check for consistency of race or gender

print(head(edstays))

write.csv(edstays, "edstays_mult_visit_with_vitals.csv")
