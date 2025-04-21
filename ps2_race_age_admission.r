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

#question: is admission affected by AGE for different race?


#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

save_filepath <- "ps2/race_age/"
dir.create(save_filepath, recursive=TRUE) 

#data on all patients
edstays <- read.csv("edstays_binary_recoded_mult_visits_mimic_with_age_vitals.csv")

#identify covariates- delete acuity, outcomes, and age
#also exclude pregnancy
vars <- colnames(edstays)
#print(vars)
covariates <- vars[!vars %in% c("is_admitted", "is_discharged",
        "X", "x", "stay_id", "time_in_ed", "age_group_50_59","age_group_70_79", "age_group_30_39", "age_group_80_89", "age_group_40_49",
        "age_group_60_69", "age_group_90", "is_asian", "is_other", "is_unknown", "is_black_african_american", 
        "is_hispanic_latino", "is_american_indian_alaska_native", "is_native_hawaiian_other_pacific_islander")]



#define age groups
ages <- c("age_group_30_39", "age_group_40_49", "age_group_50_59", "age_group_60_69", "age_group_70_79", "age_group_80_89", "age_group_90")
age_names <- c("30s", "40s", "50s", "60s", "70s", "80s", "90+")

#define races
races <- c("is_asian", "is_other", "is_unknown", "is_black_african_american", 
        "is_hispanic_latino", "is_american_indian_alaska_native", "is_native_hawaiian_other_pacific_islander", "is_white")
race_names <- c("Asian", "Other", "Unknown", "Black/African-American", "Hispanic/Latino", "American Indian", "Pacific Islander", "White")


races_to_check <- c("is_asian",  "is_black_african_american", 
        "is_hispanic_latino", "is_white")
race_names_to_check <- c("Asian", "Black/African-American", "Hispanic/Latino", "White")


#define reference group- 18-29 year olds
edstays$ref_age <- ifelse(rowSums(edstays[ages])==0, 1, 0)


#define reference group
edstays$is_white <- ifelse(rowSums(edstays[races[!(races == "is_white")]])==0, 1, 0)


#create a dataframe to save odds ratios which result from different methods
na_vector <- c(rep(NA, 2*length(ages)*length(races_to_check)))
odds_ratios <- data.frame(race=na_vector, age=na_vector, method=na_vector, estimate=na_vector, lower_conf=na_vector, upper_conf=na_vector,
        pval=na_vector, num_treated=na_vector, num_control=na_vector)
index <- 1

#threshold for caliper matching
threshold <- 0.2

#define the outcome
outcome_var <- "is_admitted"

for (race in races_to_check) {
    for (age in ages) {

        #define the problem
        adjustment_var <- age

        #do not balance on pregnancy for age 50+
        if (!((age=="age_group_30_39" | age=="age_group_40_49"))) {
            covariates <- covariates[!(covariates=="pregnant")]
        }

        #for odds ratio to make any sense, we need to limit comparison to the specified age and a reference group
        to_fit <- edstays %>% filter((.data[[age]]==1 | ref_age==1), .data[[race]]==1)

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
        ggsave(paste0(save_filepath, adjustment_var, "_", outcome_var, "_", threshold, "_balance_plot.pdf"), plot=p)


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
        mcnemar_result <- c(race, age, paste0("match_", threshold), as.numeric(test$estimate), as.numeric(test$conf.int[1]), as.numeric(test$conf.int[2]), as.numeric(test$p.val), length(treated_indices), length(control_indices))
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

        odds_ratios[index,] <- c(race, age, "cov_adjustment", exp(coeffs[2]), lower_confs[2], upper_confs[2], p_values[2], num_treated, num_control)
        index <- index + 1
    }
}


#now save both sets of covariates
write.csv(odds_ratios, paste0(save_filepath, outcome_var,"_odds_ratios.csv"))

#now plot odds ratios
method_names <- c("using McNemar's test on matched cohort", "using propensity score as an additional covariate")
methods <- c(paste0("match_", threshold), "cov_adjustment")

#convert to numeric
odds_ratios$estimate <- as.numeric(odds_ratios$estimate)
odds_ratios$lower_conf <- as.numeric(odds_ratios$lower_conf)
odds_ratios$upper_conf <- as.numeric(odds_ratios$upper_conf)


#different plot for each method
for (k in 1:length(races_to_check)) {
    for (i in 1:length(methods)) {
        ors <- odds_ratios %>% filter(method==methods[i], race==races_to_check[k])
        ors$odds_label <- paste0(
            round(ors$estimate, 2), 
            " (", 
            round(ors$lower_conf, 2), 
            "-", 
            round(ors$upper_conf, 2), 
            ")")

        ors$age <- factor(ors$age, levels=ages, labels=age_names)
        y_labels_for_full_adjustment <- paste0(age_names, "\n(N=", ors$num_treated, ")")

        p <- ggplot(ors, aes(y = age, x = estimate)) +
            geom_point(size = 3) + 
            theme(text = element_text(size =17)) +
            geom_errorbar(aes(xmin = lower_conf, xmax = upper_conf), width = 0.2) +
            geom_vline(xintercept = 1, linetype = "solid", color = "black", alpha=0.5) + 
            geom_text(aes(label = odds_label), hjust = 0.5, vjust = -1.7, size = 4) +
            labs(title=paste("Effect of age on likelihood of admitted to hospital, for people of", race_names_to_check[k], "race") , subtitle=paste0(method_names[i], ", relative to 18-29 year olds"), x = "Odds Ratios (95% CI)", y="") +
            theme_minimal() +
            scale_y_discrete(labels = y_labels_for_full_adjustment) +
            theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))


        ggsave(paste0(save_filepath, races_to_check[k], "_",  outcome_var, "_", methods[i], "_ors.pdf"), plot=p)
    }
}