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

#question: are WAIT TIMES affected by RACE+GENDER?



#here we're using binary outcomes-- does a patient have to wait more than X hours?
wait_times <- c(6.5*60, 9*60) #in minutes
wait_names <- c("top_50", "top_25")

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

save_filepath <- "ps2/rg_intersection/wait_times/admitted/"
dir.create(save_filepath, recursive=TRUE) 

#data on all patients
adm_edstays <- read.csv("adm_edstays_binary_recoded_mult_visits_mimic_with_age_vitals.csv")

#filter out those who wait more than 24 hours
adm_edstays <- adm_edstays %>% filter(time_in_ed <= 24*60)


#identify covariates- delete outcomes, race and gender
#we are comparing women of color to white men, so delete gyn. complaints
vars <- colnames(adm_edstays)
covariates <- vars[!vars %in% c("is_admitted", "is_discharged",  
        "X", "x", "stay_id", "time_in_ed", "is_asian", "is_other", "is_unknown", "is_black_african_american", 
        "is_hispanic_latino", "is_american_indian_alaska_native", "is_native_hawaiian_other_pacific_islander", "is_male",
        "pregnant", "vaginal_bleeding")]


#add columns indicating 'long wait'
for (j in 1:length(wait_times)) {

    wait_time <- wait_times[j]
    wait_name <- wait_names[j]

    adm_edstays[[wait_name]] <- ifelse(adm_edstays$time_in_ed > wait_time, 1, 0)

}

#define races
races <- c("is_asian", "is_other", "is_unknown", "is_black_african_american", 
        "is_hispanic_latino", "is_american_indian_alaska_native", "is_native_hawaiian_other_pacific_islander")
race_names <- c("Asian", "Other", "Unknown", "Black/African-American", "Hispanic/Latino", "American Indian", "Pacific Islander")



#define reference group
adm_edstays$is_white <- ifelse(rowSums(adm_edstays[races])==0, 1, 0)

#create a dataframe to save odds ratios which result from different methods
na_vector <- c(rep(NA, 2*length(races)*length(wait_times)))
odds_ratios <- data.frame(race=na_vector, wait_time=na_vector, method=na_vector, estimate=na_vector, lower_conf=na_vector, upper_conf=na_vector,
        pval=na_vector, num_treated=na_vector, num_control=na_vector)
index <- 1

#threshold for caliper matching
threshold <- 0.2


for (race in races) {

    #define the problem- this refers to WOMEN of color
    adjustment_var <- race

    #for odds ratio to make any sense, we need to limit comparison to the specified race (women) and a reference group (white men)
    to_fit <- adm_edstays %>% filter((.data[[race]]==1 & is_male==0) | (is_white==1 & is_male==1))

    #calculate propensity score and add to table
    ps_formula <- as.formula(paste(adjustment_var, "~", paste(covariates, collapse = " + ")))
    ps_model <- glm(ps_formula, data=to_fit, family="binomial"(link="logit"))
    to_fit$propensity <- ps_model$fitted.values



    #Method 1: matching with caliper threshold

    #formula is supplied to this object for balance checking purposes; distance overrides it, so PS is not calculated twice
    match_object <- matchit(ps_formula, data=to_fit, method = 'nearest', caliper=threshold*sd(to_fit$propensity), distance = to_fit$propensity)


    #print out balance tables
    p <- love.plot(match_object, vars = covariates, stats = c("mean.diffs"),
            thresholds = c(m = .1, v = 2), abs = TRUE, 
            binary = "std",
            var.order = "unadjusted",
            stars = "std"
            ) + theme(axis.text.y = element_text(size = 8))
    ggsave(paste0(save_filepath, adjustment_var, "_", threshold, "_balance_plot.pdf"), plot=p)

    for (j in 1:length(wait_names)) {

        wait_name <- wait_names[j]
        wait_time <- wait_times[j]

        outcome_var <- wait_name

        #now run McNemar's test
        matched <- as.data.table(match_data(match_object))

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


        #conduct McNemar's test
        matrix <- matrix(c(a, b, c, d), 2, 2)
        test <- mcnemar.exact(matrix)

        #save in table
        mcnemar_result <- c(race, wait_time, paste0("match_", threshold), as.numeric(test$estimate), as.numeric(test$conf.int[1]), as.numeric(test$conf.int[2]), as.numeric(test$p.val), length(treated_indices), length(control_indices))
        odds_ratios[index,] <- mcnemar_result
        index <- index + 1

        #Method 2: using propensity score as an additional covariate
        adj_formula <- as.formula(paste(outcome_var, " ~ " , adjustment_var, " + propensity + ", paste(covariates, collapse = " + ")))
        adj_model <- glm(adj_formula, family="binomial"(link="logit"), data=to_fit)

        coeffs <- adj_model$coefficients
        conf_ints <- data.frame(exp(confint.default(adj_model)))
        conf_ints <- conf_ints@.Data
        lower_confs <- conf_ints[[1]]
        upper_confs <- conf_ints[[2]]
        p_values <- coef(summary(adj_model))[,4]

        #measure treated and control here
        num_treated <- length(which(to_fit[[adjustment_var]]==1))
        num_control <- length(which(to_fit[[adjustment_var]]==0))

        odds_ratios[index,] <- c(race, wait_time, "cov_adjustment",  exp(coeffs[2]), lower_confs[2], upper_confs[2], p_values[2], num_treated, num_control)
        index <- index + 1
    }
}

#now save both sets of OBSERVED odds ratios
write.csv(odds_ratios, paste0(save_filepath, "observed_odds_ratios.csv"))

#calculate EXPECTED odds ratios (if race and gender interact independently)
gender_ors <- read.csv(paste0("ps2/gender/wait_times/admitted/wait_time_odds_ratios.csv"))
race_ors <- read.csv(paste0("ps2/race/wait_times/admitted/odds_ratios.csv"))

#define methods
method_names <- c("using McNemar's test on matched cohort", "using propensity score as an additional covariate")
methods <- c(paste0("match_", threshold), "cov_adjustment")

expected_race_gender_ors <- data.frame(race=na_vector, wait_time=na_vector, method=na_vector, estimate=na_vector, lower_conf=na_vector, upper_conf=na_vector)

index <- 1

for (i in 1:length(methods)) {
    for (k in 1:length(wait_times)) {
        g_index <- which(gender_ors$method==methods[i] & gender_ors$wait_time==wait_times[k])[1]
        male_est <- gender_ors$estimate[g_index]
        male_lc <- gender_ors$lower_conf[g_index]
        male_uc <- gender_ors$upper_conf[g_index]

        #convert to log scale and invert
        women_log_est <- -log(male_est)
        women_std <- abs(log(male_uc)-log(male_lc))/(2*1.96) #standard error

        #now calculate expected odds ratio
        for (j in 1:length(races)) {
            r_index <- which(race_ors$method==methods[i] & race_ors$race==races[j] & race_ors$wait_time==wait_times[k])[1]
            race_est <- race_ors$estimate[r_index]
            race_lc <- race_ors$lower_conf[r_index]
            race_uc <- race_ors$upper_conf[r_index]

            race_log_est <- log(race_est)
            race_std <- abs(log(race_uc)-log(race_lc))/(2*1.96)

            #combine standard errors per Cam Physics advice
            #since these are standard errors in things we are adding together (on log scale)

            combined_std <- sqrt(race_std**2 + women_std**2)

            exp_or <- exp(race_log_est+women_log_est)
            exp_lc <- exp(race_log_est+women_log_est - 1.96*combined_std)
            exp_uc <- exp(race_log_est+women_log_est + 1.96*combined_std)

            expected_race_gender_ors[index,] <- c(races[j], wait_times[k], methods[i], exp_or, exp_lc, exp_uc)
            index <- index + 1
        }
    }

}

#now save both sets of EXPECTED odds ratios
write.csv(expected_race_gender_ors, paste0(save_filepath, "expected_odds_ratios.csv"))

for (i in 1:length(methods)) {
    for (j in 1:length(wait_names)) {
        observed_ors <- odds_ratios %>% filter(method==methods[i], wait_time==wait_times[j])
        expected_ors <- expected_race_gender_ors %>% filter(method==methods[i], wait_time==wait_times[j])


        y_labels_for_full_adjustment <- paste0(race_names, "\n(N=", observed_ors$num_treated, ")")
        observed_ors <- subset(observed_ors, select=c("race", "wait_time", "method", "estimate", "lower_conf", "upper_conf"))

        print(head(expected_ors))
        print(head(observed_ors))

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
        labs(title=paste("Odds of waiting longer than", wait_times[j], "minutes to be admitted, for women of color"), subtitle=paste0(method_names[i], ", relative to white men"), x = "Odds Ratios (95% CI)", y="") +
        theme_minimal() +
        scale_color_manual(values = c("Expected" = "blue", "Observed" = "black")) + 
        theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

        ggsave(paste0(save_filepath, wait_times[j], "_", methods[i], "_interaction_plot.pdf"), plot=p)

    }
}