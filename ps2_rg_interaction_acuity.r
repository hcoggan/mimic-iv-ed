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
library(MatchIt)
library(data.table)
library(exact2x2)
library(cobalt)


#second attempt at propensity scoring
#two methods: propensity matching with a caliper distance of 0.2*standard deviation
#and adjustment with the propensity score as an additional covariate

#question: is URGENCY affected by RACE+GENDER?
#here we calculate expected odds ratios and make plots


#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

save_filepath <- "ps2/rg_intersection/"
dir.create(save_filepath, recursive=TRUE) 

#data on all patients
edstays <- read.csv("edstays_binary_recoded_mult_visits_mimic_with_age.csv")

#identify covariates- delete acuity, outcomes, race and gender
#we are comparing women of color to white men, so delete gyn. complaints
vars <- colnames(edstays)
covariates <- vars[!vars %in% c("is_admitted", "is_discharged", "acuity_high", "acuity_very_high", "acuity_low", "acuity_unassigned", 
        "X", "x", "stay_id", "time_in_ed", "is_asian", "is_other", "is_unknown", "is_black_african_american", 
        "is_hispanic_latino", "is_american_indian_alaska_native", "is_native_hawaiian_other_pacific_islander", "is_male",
        "pregnant", "vaginal_bleeding")]

#recode previous visits and previous admissions as binary
edstays$previous_visits <- ifelse(edstays$previous_visits>0, 1, 0)
edstays$previous_admissions <- ifelse(edstays$previous_admissions>0, 1, 0)

#create a combined column for ESI scores 1/2- this is our outcome
edstays$urgent <- edstays$acuity_high + edstays$acuity_very_high

#define races
races <- c("is_asian", "is_other", "is_unknown", "is_black_african_american", 
        "is_hispanic_latino", "is_american_indian_alaska_native", "is_native_hawaiian_other_pacific_islander")
race_names <- c("Asian", "Other", "Unknown", "Black/African-American", "Hispanic/Latino", "American Indian", "Pacific Islander")


#define reference group
edstays$is_white <- ifelse(rowSums(edstays[races]==0), 1, 0)


#calculate EXPECTED odds ratios (if race and gender interact independently)
gender_ors <- read.csv(paste0("ps2/gender/",outcome_var,"_odds_ratios.csv"))
race_ors <- read.csv(paste0("ps2/race/",outcome_var,"_odds_ratios.csv"))

#define methods
method_names <- c("using McNemar's test on matched cohort", "using propensity score as an additional covariate")
methods <- c(paste0("match_", threshold), "cov_adjustment")

expected_race_gender_ors <- data.frame(race=na_vector, method=na_vector, estimate=na_vector, lower_conf=na_vector, upper_conf=na_vector)

index <- 1
for (i in 1:length(methods)) {
    g_index <- which(gender_ors$method==methods[i])[1]
    male_est <- gender_ors$estimate[g_index]
    male_lc <- gender_ors$lower_conf[g_index]
    male_uc <- gender_ors$upper_conf[g_index]

    #convert to log scale and invert
    women_log_est <- -log(male_est)
    women_std <- (log(male_uc)-log(male_lc))/2 #this is 95%-erf multiplied by standard error in gender estimation

    #now calculate expected odds ratio
    for (j in 1:length(races)) {
        r_index <- which(race_ors$method==methods[i] & race_ors$race==races[j])[1]
        race_est <- race_ors$estimate[r_index]
        race_lc <- race_ors$lower_conf[r_index]
        race_uc <- race_ors$upper_conf[r_index]

        race_log_est <- log(race_est)
        race_std <- (log(race_uc)-log(race_lc))/2

        exp_or <- exp(race_log_est+women_log_est)
        exp_lc <- exp(race_log_est+women_log_est - (race_std + women_std))
        exp_uc <- exp(race_log_est+women_log_est + (race_std + women_std))

        expected_race_gender_ors[index,] <- c(races[j], methods[i], exp_or, exp_lc, exp_uc)
        index <- index + 1
    }

}
#now save both sets of EXPECTED odds ratios
write.csv(expected_race_gender_ors, paste0(save_filepath, outcome_var,"_expected_odds_ratios.csv"))


for (i in 1:length(methods)) {
    observed_ors <- odds_ratios %>% filter(method==methods[i])
    expected_ors <- expected_race_gender_ors %>% filter(method==methods[i])
    observed_ors <- subset(observed_ors, select=c("race", "method", "estimate", "lower_conf", "upper_conf"))
    expected_ors$Ratio <- "Expected"
    observed_ors$Ratio <- "Observed"
    combined_data <- rbind(expected_ors, observed_ors)

    combined_data$estimate <- as.numeric(combined_data$estimate)
    combined_data$lower_conf <- as.numeric(combined_data$lower_conf)
    combined_data$upper_conf <- as.numeric(combined_data$upper_conf)


    #now make a plot
    combined_data$race <- factor(combined_data$race, levels=races, labels=race_names)

    p <- ggplot(combined_data, aes(y = race, x = estimate, color = Ratio)) +
        geom_point(size = 3) + 
        theme(text = element_text(size =17)) +
        geom_errorbar(aes(xmin = lower_conf, xmax = upper_conf), width = 0.2) +
        geom_vline(xintercept = 1, linetype = "solid", color = "black", alpha=0.5) + 
        labs(title="Odds of being triaged as 'urgent' for women of color", subtitle=paste(method_names[i], ", relative to white men"), x = "Odds Ratios (95% CI)", y="") +
        theme_minimal() +
        scale_color_manual(values = c("Expected" = "blue", "Observed" = "black")) + 
        theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

    ggsave(paste0(save_filepath, outcome_var, "_", methods[i], "_interaction_plot.pdf"), plot=p)

}