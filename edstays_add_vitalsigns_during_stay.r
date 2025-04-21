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

#File adds vital signs collected during stay, with a flag for whether or not a category was ever met

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

#edstays has already undergone binary recode
edstays <- read.csv("edstays_binary_recoded_mult_visits_mimic_with_age.csv")
adm_edstays <- read.csv("adm_edstays_binary_recoded_mult_visits_mimic_with_age.csv")
vitalsigns <-  read.csv("mimic-iv-ed-2.2/ed/vitalsign.csv")


print("dataframes loaded")


#Apply same filtration to vitalsigns as we do to edstays


#first, if something is clearly intended as a number between 0 and 10, keep it
valid_pain_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")


#finally, condense every integer up to 20 down to 10
indices <- which((!is.na(vitalsigns$pain)) & (as.numeric(vitalsigns$pain) < 21) & (as.numeric(vitalsigns$pain) > 10))
vitalsigns$pain <- replace(vitalsigns$pain, indices, "10")

for (pain in valid_pain_levels) {
    #find those with random punctuation either side
    intended_level <- which((vitalsigns$pain == paste(" ", pain, sep="")) |
         (vitalsigns$pain == paste(pain, " ", sep="")) | (vitalsigns$pain == paste("  ", pain, sep="")) | 
         (vitalsigns$pain == paste(pain, "  ", sep="")) | (vitalsigns$pain == paste(pain, ".", sep="")) | 
         (vitalsigns$pain == paste(".", pain, sep="")) | (vitalsigns$pain == paste(pain, "-", sep="")) | 
         (vitalsigns$pain == paste("-", pain, sep="")) | (vitalsigns$pain == paste(pain, "/10", sep="")) |
         (vitalsigns$pain == paste(pain, ",", sep="")) | (vitalsigns$pain == paste("0", pain, sep="")) | (vitalsigns$pain == paste(",", pain, sep="")))
    vitalsigns$pain <- replace(vitalsigns$pain, intended_level, pain)
}

#first condense intermediate or out of range values
for (pain in valid_pain_levels) {
    for (pain_ in valid_pain_levels) {
        p <- as.numeric(pain)
        p_ <- as.numeric(pain_)
        #find ranges (i.e. '6-7')
        rep_value <- (p+p_)/2
        low_int <- floor(rep_value) #lowest integer
        indices <- which(vitalsigns$pain == paste(pain, pain_, sep="-"))
        if (!(rep_value == low_int)) { #if the middle of this range is not an integer
            prob_higher_integer <- rep_value - low_int #i.e. 3.8 is assigned to 4 80% of the time, and to 3 20% of the time
            for (index in indices) {
                if (runif(1) < prob_higher_integer) { #assign the lower or higher integer probabilistically 
                    vitalsigns$pain[index] <- as.character(low_int + 1)
                } else {
                    vitalsigns$pain[index] <- as.character(low_int)
                }
            }
        } else { #middle of this range is an integer
                for (index in indices) {
                    vitalsigns$pain[index] <- as.character(low_int)
                }
        }
        #now find values between 20 and 100
        val <- 10*p+p_
        if (val > 20) {
            val_char <- as.character(val)
            indices <- which(vitalsigns$pain == val_char)
            #replace it assuming it was meant to be between 1 and 100
            rep_value <- val/10 
            low_int <- floor(rep_value) #lowest integer
            if (!(rep_value == low_int)) { #if the middle of this range is not an integer
                prob_higher_integer <- rep_value - low_int #i.e. 3.8 is assigned to 4 80% of the time, and to 3 20% of the time
                for (index in indices) {
                    if (runif(1) < prob_higher_integer) { #assign the lower or higher integer probabilistically 
                        vitalsigns$pain[index] <- as.character(low_int + 1)
                    } else {
                        vitalsigns$pain[index] <- as.character(low_int)
                    }
                }
            } else { #middle of this range is an integer
                for (index in indices) {
                    vitalsigns$pain[index] <- as.character(low_int)
                }
        }
        }
        #now find intermediate values e.g. 6.5
        val_char <- paste(pain, pain_, sep=".")
        rep_value <- as.numeric(val_char)
        low_int <- floor(rep_value) #integer below it
        indices <- which(vitalsigns$pain == val_char)
        prob_higher_integer <- rep_value - low_int #i.e. 3.8 is assigned to 4 80% of the time, and to 3 20% of the time
        for (index in indices) {
            if (runif(1) < prob_higher_integer) { #assign the lower or higher integer probabilistically 
                vitalsigns$pain[index] <- as.character(low_int + 1)
            } else {
                vitalsigns$pain[index] <- as.character(low_int)
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

unable_indices <- which(vitalsigns$pain %in% unable_synonyms)
vitalsigns$pain <- replace(vitalsigns$pain, unable_indices, "UTA")
                       



#final values
valid_pain_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "UTA")

#condense everything else down to 'Other'
other_indices <- which(!(vitalsigns$pain %in% valid_pain_levels))
vitalsigns$pain <- replace(vitalsigns$pain, other_indices, "non-numeric")

#now condense numeric response into none, mild, moderate, severe
none <- which(vitalsigns$pain=="0")
mild <- which((vitalsigns$pain=="1") | (vitalsigns$pain=="2") | (vitalsigns$pain=="3"))
moderate <- which((vitalsigns$pain=="4") | (vitalsigns$pain=="5") | (vitalsigns$pain=="6") | (vitalsigns$pain=="7"))
severe <- which((vitalsigns$pain=="8") | (vitalsigns$pain=="9") | (vitalsigns$pain=="10"))

vitalsigns$pain <- replace(vitalsigns$pain, none, "none")
vitalsigns$pain <- replace(vitalsigns$pain, mild, "mild")
vitalsigns$pain <- replace(vitalsigns$pain, moderate, "moderate")
vitalsigns$pain <- replace(vitalsigns$pain, severe, "severe")

#filter temperature, resp rate, blood pressure on La Cava Github


vitalsigns <- vitalsigns %>% filter(between(temperature, 95, 105) & between(resprate, 2, 200)
        & between(o2sat, 50, 100) & between(sbp, 30, 400) & between(dbp,30, 300)
     & between(heartrate, 30, 300))

#define fever per Andrew Fine and Pradip Chaudhari
fever_grade <- c(rep("none", nrow(vitalsigns)))
fever <- which(vitalsigns$temperature > 100.4)
not_known <- which(is.na(vitalsigns$temperature))
fever_grade[fever] <- "fever"
fever_grade[not_known] <- "not_known"

vitalsigns$fever_grade <- fever_grade

#define O2 level per Wold 2015 (and Jaya Aysola)
o2_level <- c(rep("normal", nrow(vitalsigns)))
low <- which(vitalsigns$o2sat <= 95.0 & vitalsigns$o2sat > 85)
very_low <- which(vitalsigns$o2sat <= 85)
not_known <- which(is.na(vitalsigns$o2sat))
#very_low <- which(vitalsigns$o2sat <= 92)


o2_level <- replace(o2_level, low, "85_to_95")
o2_level <- replace(o2_level, very_low, "below_85")
o2_level <- replace(o2_level, not_known, "not_known")
#o2_level <- replace(o2_level, very_low, "very low")

vitalsigns$o2_level <- o2_level


#define resp level per NIH
resp_level <- c(rep("normal", nrow(vitalsigns)))
low <- which(vitalsigns$resprate < 12)
high <- which(vitalsigns$resprate > 20)
not_known <- which(is.na(vitalsigns$resp_level))


resp_level <- replace(resp_level, low, "low")
resp_level <- replace(resp_level, high, "high")
resp_level <- replace(resp_level, not_known, "not_known")


vitalsigns$resp_level <- resp_level

#also heart rate levels, per Nursing Skills textbook
heart_level <- c(rep("normal", nrow(vitalsigns)))
low <- which(vitalsigns$heartrate < 60)
high <- which(vitalsigns$heartrate > 100 & vitalsigns$heartrate <= 150)
very_high <- which(vitalsigns$heartrate > 150)
not_known <- which(is.na(vitalsigns$heartrate))


heart_level <- replace(heart_level, low, "low")
heart_level <- replace(heart_level, high, "high")
heart_level <- replace(heart_level, very_high, "very_high")
heart_level <- replace(heart_level, not_known, "not_known")


vitalsigns$heart_level <- heart_level

#also blood pressure levels, per Johns Hopkins

bp_level <- c(rep("normal", nrow(vitalsigns)))
low <- which((vitalsigns$sbp < 90) & (vitalsigns$dbp < 60))
elevated <- which((vitalsigns$sbp >= 120) & (vitalsigns$sbp < 130) & (vitalsigns$dbp < 80))
stage_1_high <- which(((vitalsigns$sbp >= 130) & (vitalsigns$sbp < 140)) | ((vitalsigns$dbp >= 80) & (vitalsigns$dbp < 90)))
stage_2_high <- which(((vitalsigns$sbp >= 140) & (vitalsigns$sbp < 180)) | (vitalsigns$dbp >= 90))
unstable <- which(vitalsigns$sbp >= 180)
not_known <- which(is.na(vitalsigns$sbp) | is.na(vitalsigns$dbp))

bp_level <- replace(bp_level, low, "low")
bp_level <- replace(bp_level, elevated, "elevated")
bp_level <- replace(bp_level, stage_1_high, "stage_1_high")
bp_level <- replace(bp_level, stage_2_high, "stage_2_high") #140/90 is stage 2 hypertension threshold-- if something 'could also be' stage 1 by above definition, we assign it stage 2
bp_level <- replace(bp_level, unstable, "unstable") 
bp_level <- replace(bp_level, not_known, "not_known")



#check that everyone not caught by these does in fact have normal blood pressure 
classified_normal <- which(bp_level == "normal")
actually_normal <- which((vitalsigns$sbp < 120) & (vitalsigns$dbp < 80))
print("Set difference") #we want this to be 0
print(length(setdiff(classified_normal, actually_normal)))

vitalsigns$bp_level <- bp_level
#now: for each nonreference category, ask whether it is EVER true of a patient through their stay in the ED
#exclude 'not known'
#eg is a fever EVER recorded?

coded_data <- data.frame(stay_id=vitalsigns$stay_id)

vital_names <- c("pain", "fever_grade", "o2_level", "resp_level", "heart_level", "bp_level")


for (name in vital_names) {
    values <- unique(vitalsigns[[name]]) #values of the vital category
    values <- values[!(values %in% c("not_known", "none", "normal"))] #exclude 'not recorded', for all except pain
    for (value in values) {
        col_name <- paste0("during_stay_", name, "_", value)
        coded_data[[col_name]] <- ifelse(vitalsigns[[name]]==value, 1, 0)
    }
}

print("success")
during_stay_vars <- colnames(coded_data)
during_stay_vars <- during_stay_vars[!(during_stay_vars == "stay_id")]
print(during_stay_vars)

#create an empty dataframe
flag_data <- data.frame(stay_id=unique(coded_data$stay_id))
for (var in during_stay_vars) {
    flag_data[[var]] <- NA
}

num <- length(during_stay_vars)

#takes an ID and flags whether any 1s are witnessed
check_vitals_during_stay <- function(id) {
    indices <- which(coded_data$stay_id==id)
    vec <- c(rep(NA, num))
    for (j in 1:length(during_stay_vars)) {
        vec[j] <- ifelse(any(coded_data[indices, during_stay_vars[j]]==1), 1, 0)
    }
    return(vec)
}

ids <- unique(coded_data$stay_id)

flag_data <- purrr::map(ids, check_vitals_during_stay, .progress=TRUE)
flag_data <- as.data.frame(do.call(rbind, flag_data))


colnames(flag_data) <- c(during_stay_vars)
#rownames(flag_data) <- seq(1:nrow(flag_data))

flag_data$stay_id <- ids



edstays <- edstays %>% left_join(flag_data, by="stay_id", copy=TRUE)
adm_edstays <- adm_edstays %>% left_join(flag_data, by="stay_id", copy+TRUE)

#if no vitals are recorded, just zero out the flags
for (name in during_stay_vars) {
    indexes <- is.na(edstays[[name]])
    edstays[indexes, name] <- 0
    indexes <- is.na(adm_edstays[[name]])
    adm_edstays[indexes, name] <- 0
}

write.csv(edstays, "edstays_binary_recoded_mult_visits_mimic_with_age_vitals.csv")

write.csv(adm_edstays, "adm_edstays_binary_recoded_mult_visits_mimic_with_age_vitals.csv")
