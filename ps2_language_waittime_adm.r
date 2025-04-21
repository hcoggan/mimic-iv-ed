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

#Here we're trying to use propensity scoring to identify the effect of primary_language being not english on wait times
#for admitted patients

#here we're using binary outcomes-- does a patient have to wait more than X hours?
wait_times <- c(6.5*60, 9*60) #in minutes
wait_names <- c("top_50", "top_25")

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

adm_edstays <- read.csv("adm_edstays_binary_recoded_mult_visits_mimic_with_age_vitals.csv")

save_filepath <- "ps2/language/wait_times/admitted/"
dir.create(save_filepath, recursive=TRUE) 


#print(head(adm_edstays))

#sanity check: these are admitted patients, right?
assert(all(adm_edstays$is_admitted==1))

#filter out those who wait more than 24 hours
adm_edstays <- adm_edstays %>% filter(time_in_ed <= 24*60)

#identify covariates- delete acuity, outcomes, and gender
#also delete gyn. complaints
vars <- colnames(adm_edstays)
covariates <- vars[!vars %in% c("is_admitted", "is_discharged", 
        "X", "x", "stay_id", "time_in_ed", "primary_english")]

#print(covariates)



#define the problem
adm_edstays$non_english <- ifelse(adm_edstays$primary_english==0, 1, 0)
adjustment_var <- "non_english"

#add columns indicating 'long wait'
for (j in 1:length(wait_times)) {

    wait_time <- wait_times[j]
    wait_name <- wait_names[j]

    adm_edstays[[wait_name]] <- ifelse(adm_edstays$time_in_ed > wait_time, 1, 0)

}



#create a dataframe to save odds ratios which result from different methods AND WAIT TIMES
na_vector <- c(rep(NA, 2*length(wait_times)))
odds_ratios <- data.frame(method=na_vector, wait_times=na_vector, estimate=na_vector, lower_conf=na_vector, upper_conf=na_vector,
        pval=na_vector, num_treated=na_vector, num_control=na_vector)
index <- 1



#calculate propensity score and add to table
ps_formula <- as.formula(paste(adjustment_var, "~", paste(covariates, collapse = " + ")))
ps_model <- glm(ps_formula, data=adm_edstays, family="binomial"(link="logit"))
adm_edstays$propensity <- ps_model$fitted.values

#Method 1: matching with caliper threshold
threshold <- 0.2


#formula is supplied to this object for balance checking purposes; distance overrides it, so PS is not calculated twice
match_object <- matchit(ps_formula, data=adm_edstays, method = 'nearest', caliper=threshold*sd(adm_edstays$propensity), distance = adm_edstays$propensity)


#print out balance tables
p <- love.plot(match_object, vars = covariates, stats = c("mean.diffs"),
          thresholds = c(m = .1, v = 2), abs = TRUE, 
          binary = "std",
          var.order = "unadjusted",
          stars = "std"
          ) + theme(axis.text.y = element_text(size = 8))
ggsave(paste0(save_filepath, "_", threshold, "_balance_plot.pdf"), plot=p)

#now run McNemar's test on these three outcomes
matched <- as.data.table(match_data(match_object))

#iterate over the wait times
for (j in 1:length(wait_times)) {

    outcome_var <- wait_names[j]
    wait_time <- wait_times[j]

    #make two data tables for treated and untreated data, rejoin on pair
    treated_indices <- which(matched[[adjustment_var]]==1)
    control_indices <- which(matched[[adjustment_var]]==0)
    treated_df <- data.frame(pair_index=matched[treated_indices, subclass], treatment_outcome=matched[[outcome_var]][treated_indices])
    control_df <- data.frame(pair_index=matched[control_indices, subclass], control_outcome=matched[[outcome_var]][control_indices])
    overall_table <- inner_join(treated_df, control_df, by="pair_index") #

    a <- sum(overall_table$treatment_outcome==1 & overall_table$control_outcome==1)
    c <- sum(overall_table$treatment_outcome==1 & overall_table$control_outcome==0)
    b <- sum(overall_table$treatment_outcome==0 & overall_table$control_outcome==1)
    d <- sum(overall_table$treatment_outcome==0 & overall_table$control_outcome==0)

    #now print: no pictures here
    matrix <- matrix(c(a, b, c, d), 2, 2)
    test <- mcnemar.exact(matrix)

    #save in table
    mcnemar_result <- c(paste0("match_", threshold), wait_time, test$estimate, test$conf.int[1], test$conf.int[2], test$p.val, length(treated_indices), length(control_indices))
    odds_ratios[index,] <- mcnemar_result
    index <- index + 1

    #Method 2: using propensity score as an additional covariate
    adj_formula <- as.formula(paste(outcome_var, " ~ " , adjustment_var, " + propensity + ", paste(covariates, collapse = " + ")))
    adj_model <- glm(adj_formula, family="binomial"(link="logit"), data=adm_edstays)

    coeffs <- adj_model$coefficients
    conf_ints <- data.frame(exp(confint.default(adj_model)))
    conf_ints <- conf_ints@.Data
    lower_confs <- conf_ints[[1]]
    upper_confs <- conf_ints[[2]]
    p_values <- coef(summary(adj_model))[,4]

    #measure treated and control here
    num_treated <- length(which(adm_edstays[[adjustment_var]]==1))
    num_control <- length(which(adm_edstays[[adjustment_var]]==0))

    odds_ratios[index,] <- c("cov_adjustment", wait_time, exp(coeffs[2]), lower_confs[2], upper_confs[2], p_values[2], num_treated, num_control)
    index <- index + 1
}

#now save both sets of covariates
write.csv(odds_ratios, paste0(save_filepath,"wait_time_odds_ratios.csv"))
