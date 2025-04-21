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

#This file recodes all exposures  as binary variables; only acts on adm_edstays, filtered by first visit with date of death

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

adm_edstays <- read.csv("adm_edstays_first_visit_with_linked_death.csv")

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
        varname <- paste("fever_", level)
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

#primary-english is already binary
coded_data$primary_english <- adm_edstays$primary_english

#we don't need admission and discharged; they were all admitted

#previous visits are non-binary variables-- we can look at impact of all previous visits

coded_data$previous_visits <- adm_edstays$previous_visits
coded_data$previous_admissions <- adm_edstays$previous_admissions

#copy across complaint variables
#load complaints alone
complaints <- read.csv("binarised_complaints_alone.csv")

for (name in colnames(complaints)) {
    if (name != "stay_id") {
        coded_data[[name]] <- adm_edstays[[name]]
    }
}

coded_data$time_in_ed <- adm_edstays$time_in_ed
coded_data$time_to_death <- adm_edstays$time_to_death

#keep stay_id
coded_data$stay_id <- adm_edstays$stay_id


#clean up names
coded_data <- clean_names(coded_data)

print(head(coded_data))

write.csv(coded_data, "adm_edstays_first_visit_with_linked_death_binary_recode.csv")

