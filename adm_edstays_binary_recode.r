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
hgd()

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(janitor)

#This file recodes all exposures  as binary variables; only acts on adm_edstays (includes extra data for admitted patients) is dealt with in a different file

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

adm_edstays <- read.csv("adm_edstays.csv")

print(head(adm_edstays))




#create a dataframe
num_patients <- nrow(adm_edstays)

is_male_recode <- c("M"=1, "F"=0)
is_male <- dplyr::recode(adm_edstays$gender, !!! is_male_recode)

coded_data <- data.frame(is_male=is_male) #0 corresponds to reference (female)



#turn each race into a binary is_Black, etc. 0 corresponds to reference (white)
ref_race <- "White" 
races <- unique(adm_edstays$race)
for (race in races) {
    if (!(race==ref_race)) {
        varname <- paste0("is_", race)
        values <- c(rep(0, num_patients))
        indices <- which(adm_edstays$race == race)
        values <- replace(values, indices, 1)
        coded_data[[varname]] <- values
    }
}

#same for modes
ref_mode <- "WALK IN"
modes <- unique(adm_edstays$arrival_transport)
for (mode in modes) {
    if (!(mode == ref_mode)) {
        varname <- mode
        values <- c(rep(0, num_patients))
        indices <- which(adm_edstays$arrival_transport == mode)
        values <- replace(values, indices, 1)
        coded_data[[varname]] <- values
    }
}

#now acuity
#reference level is moderate; create binary variables corresponding to all other levels
ref_level <- "moderate"
for (level in unique(adm_edstays$acuity)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$acuity==level)
        varname <- paste("acuity_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#similarly, treat all pain as its own category
ref_level <- "none"
for (level in unique(adm_edstays$pain)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$pain==level)
        varname <- paste("pain_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now look at fever
ref_level <- "none"
for (level in unique(adm_edstays$fever_grade)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$fever_grade==level)
        varname <- paste("temp_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}


#now O2 level
ref_level <- "normal"
for (level in unique(adm_edstays$o2_level)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$o2_level==level)
        varname <- paste("o2_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now O2 level
ref_level <- "normal"
for (level in unique(adm_edstays$o2_level)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$o2_level==level)
        varname <- paste("o2_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}
#now resp level
ref_level <- "normal"
for (level in unique(adm_edstays$resp_level)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$resp_level==level)
        varname <- paste("resp_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now heart level
ref_level <- "normal"
for (level in unique(adm_edstays$heart_level)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$heart_level==level)
        varname <- paste("heart_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}
#now BP level
ref_level <- "normal"
for (level in unique(adm_edstays$bp_level)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$bp_level==level)
        varname <- paste("bp_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now time of day
ref_level <- "afternoon"
for (level in unique(adm_edstays$time_of_day)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$time_of_day==level)
        varname <- paste("time_of_day_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}


#now age
ref_level <- "18-29"
for (level in unique(adm_edstays$age_group)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$age_group==level)
        varname <- paste("age_group_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}


#now year group
ref_level <- "2008 - 2010"
for (level in unique(adm_edstays$anchor_year_group)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$anchor_year_group==level)
        varname <- paste("year_group_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now insurance_status
ref_level <- "Medicare"
for (level in unique(adm_edstays$insurance_group)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$insurance_group==level)
        varname <- paste("insurance_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}
#now marital status
ref_level <- "MARRIED"
for (level in unique(adm_edstays$marital_status)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$marital_status==level)
        varname <- paste("marital_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now bounceback visits
ref_level <- "no_bounceback"
for (level in unique(adm_edstays$bounceback)) {
    if (!(level==ref_level)) {
        indices <- which(adm_edstays$bounceback==level)
        varname <- paste(level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}


#primary-english is already binary
coded_data$primary_english <- adm_edstays$primary_english

#we don't need admission and discharged; they were all admitted


#code previous visits/admissions as binary

coded_data$recent_visit_without_admit <- ifelse(adm_edstays$recency == "recent_visit_without_admission", 1, 0)
coded_data$recent_admit <- ifelse(adm_edstays$recency == "recent_admission", 1, 0)

#no NAs in Charlson score
coded_data$charlson_1 <- ifelse(adm_edstays$charlson_score==1, 1, 0)
coded_data$charlson_2_or_more <- ifelse(adm_edstays$charlson_score>1, 1, 0)

#copy across complaint variables
#load complaints alone
complaints <- read.csv("binarised_complaints_alone.csv")

for (name in colnames(complaints)) {
    if (name != "stay_id") {
        coded_data[[name]] <- adm_edstays[[name]]
    }
}

# #finally ICD codes
# unique_codes <- unique(adm_edstays$icd_with_version)
# print(length(unique_codes)) #about 5k- flag the 22 most common, with over 1k visits

# tab <- table(adm_edstays$icd_with_version) %>% as.data.frame() %>% arrange(desc(Freq))

# write.csv(tab, "adm_edstays_icd_codes_with_version_by_frequency.csv")

# codes_to_flag <- tab$Var1[1:22]

# for (code in codes_to_flag) {
#     #print(head(edstays$icd_title==code))
#     coded_data[[code]] <- ifelse(adm_edstays$icd_with_version==code, 1, 0)
# }


coded_data$time_in_ed <- adm_edstays$time_in_ed
coded_data$time_in_hospital <- adm_edstays$time_in_hospital
coded_data$death_in_hospital <- adm_edstays$death_in_hospital

#keep stay_id
coded_data$stay_id <- adm_edstays$stay_id

#clean up names
coded_data <- clean_names(coded_data)

print(head(coded_data))

write.csv(coded_data, "adm_edstays_binary_recoded_mult_visits_mimic_with_age.csv")

