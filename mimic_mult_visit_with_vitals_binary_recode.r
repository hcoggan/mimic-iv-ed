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

#This file recodes all exposures (except those in the mimic data) as binary variables; only acts on edstays (i.e. for risk of acuity/admission)
#adm_edstays (includes extra data for admitted patients) is dealt with in a different file

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

edstays <- read.csv("edstays_with_binarised_complaints_and_age.csv")

print(head(edstays))




#create a dataframe
num_patients <- nrow(edstays)

is_male_recode <- c("M"=1, "F"=0)
is_male <- dplyr::recode(edstays$gender, !!! is_male_recode)

coded_data <- data.frame(is_male=is_male) #0 corresponds to reference (female)



#turn each race into a binary is_Black, etc. 0 corresponds to reference (white)
ref_race <- "White" 
races <- unique(edstays$race)
for (race in races) {
    if (!(race==ref_race)) {
        varname <- paste0("is_", race)
        values <- c(rep(0, num_patients))
        indices <- which(edstays$race == race)
        values <- replace(values, indices, 1)
        coded_data[[varname]] <- values
    }
}

#same for modes
ref_mode <- "WALK IN"
modes <- unique(edstays$arrival_transport)
for (mode in modes) {
    if (!(mode == ref_mode)) {
        varname <- mode
        values <- c(rep(0, num_patients))
        indices <- which(edstays$arrival_transport == mode)
        values <- replace(values, indices, 1)
        coded_data[[varname]] <- values
    }
}

#now acuity
#reference level is moderate; create binary variables corresponding to all other levels
ref_level <- "moderate"
for (level in unique(edstays$acuity)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$acuity==level)
        varname <- paste("acuity_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#similarly, treat all pain as its own category
ref_level <- "none"
for (level in unique(edstays$pain)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$pain==level)
        varname <- paste("pain_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now look at fever
ref_level <- "none"
for (level in unique(edstays$fever_grade)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$fever_grade==level)
        varname <- paste("temp_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now O2 level
ref_level <- "normal"
for (level in unique(edstays$o2_level)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$o2_level==level)
        varname <- paste("o2_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}
#now resp level
ref_level <- "normal"
for (level in unique(edstays$resp_level)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$resp_level==level)
        varname <- paste("resp_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now heart level
ref_level <- "normal"
for (level in unique(edstays$heart_level)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$heart_level==level)
        varname <- paste("heart_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}
#now BP level
ref_level <- "normal"
for (level in unique(edstays$bp_level)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$bp_level==level)
        varname <- paste("bp_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now age
ref_level <- "18-29"
for (level in unique(edstays$age_group)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$age_group==level)
        varname <- paste("age_group_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}


#now year group
ref_level <- "2008 - 2010"
for (level in unique(edstays$anchor_year_group)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$anchor_year_group==level)
        varname <- paste("year_group_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}

#now time of day
ref_level <- "afternoon"
for (level in unique(edstays$time_of_day)) {
    if (!(level==ref_level)) {
        indices <- which(edstays$time_of_day==level)
        varname <- paste("time_of_day_", level)
        binary <- c(rep(0, num_patients))
        binary <- replace(binary, indices, 1)
        coded_data[[varname]] <- binary
    }
}


#code previous visits/admissions as binary

coded_data$recent_visit_without_admit <- ifelse(edstays$recency == "recent_visit_without_admission", 1, 0)
coded_data$recent_admit <- ifelse(edstays$recency == "recent_admission", 1, 0)


#no NAs in Charlson score
coded_data$charlson_1 <- ifelse(edstays$charlson_score==1, 1, 0)
coded_data$charlson_2_or_more <- ifelse(edstays$charlson_score>1, 1, 0)

#copy across complaint variables
#load complaints alone
complaints <- read.csv("binarised_complaints_alone.csv")

for (name in colnames(complaints)) {
    if (name != "stay_id") {
        coded_data[[name]] <- edstays[[name]]
    }
}

# #finally ICD codes
# unique_codes <- unique(edstays$icd_with_version)
# print(length(unique_codes)) #more than 12000- flag the top 68, which have >1k visits

# tab <- table(edstays$icd_with_version) %>% as.data.frame() %>% arrange(desc(Freq))

# write.csv(tab, "icd_codes_with_version_by_frequency.csv")



# codes_to_flag <- tab$Var1[1:68]

# for (code in codes_to_flag) {
#     #print(head(edstays$icd_title==code))
#     coded_data[[code]] <- ifelse(edstays$icd_with_version==code, 1, 0)
# }


#copy outcome variables
coded_data$is_admitted <- edstays$admitted
coded_data$is_discharged <- edstays$discharged

coded_data$time_in_ed <- edstays$time_in_ed


#keep stay_id (and subject ID)
coded_data$stay_id <- edstays$stay_id

#clean up names
coded_data <- clean_names(coded_data)

print(head(coded_data))

write.csv(coded_data, "edstays_binary_recoded_mult_visits_mimic_with_age.csv")

