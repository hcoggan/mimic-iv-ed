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

#question: is ADMISSION affected by GENDER?

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

save_filepath <- "ps2/gender/"
dir.create(save_filepath, recursive=TRUE) 

#data on all patients
edstays <- read.csv("edstays_binary_recoded_mult_visits_mimic_with_age_vitals.csv")

#identify covariates- delete outcomes, and gender
#also delete gyn. complaints
vars <- colnames(edstays)
covariates <- vars[!vars %in% c("is_admitted", "is_discharged", "acuity_high",
        "X", "x", "stay_id", "time_in_ed", "pregnant", "vaginal_bleeding", "is_male")]


#define the problem
adjustment_var <- "is_male"
outcome_var <- "is_admitted"

#create a dataframe to save odds ratios which result from different methods
na_vector <- c(rep(NA, 2))
odds_ratios <- data.frame(method=na_vector, estimate=na_vector, lower_conf=na_vector, upper_conf=na_vector,
        pval=na_vector, num_treated=na_vector, num_control=na_vector)
index <- 1


#calculate propensity score and add to table
ps_formula <- as.formula(paste(adjustment_var, "~", paste(covariates, collapse = " + ")))
ps_model <- glm(ps_formula, data=edstays, family="binomial"(link="logit"))
edstays$propensity <- ps_model$fitted.values



#Method 1: matching with caliper threshold
threshold <- 0.2

#formula is supplied to this object for balance checking purposes; distance overrides it, so PS is not calculated twice
match_object <- matchit(ps_formula, data=edstays, method = 'nearest', caliper=threshold*sd(edstays$propensity), distance = edstays$propensity)


#print out balance tables
p <- love.plot(match_object, vars = covariates, stats = c("mean.diffs"),
          thresholds = c(m = .1, v = 2), abs = TRUE, 
          binary = "std",
          var.order = "unadjusted",
          stars = "std"
          ) + theme(axis.text.y = element_text(size = 8))
ggsave(paste0(save_filepath, outcome_var, "_", threshold, "_balance_plot.pdf"), plot=p)


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


#now print: no pictures here
matrix <- matrix(c(a, b, c, d), 2, 2)
test <- mcnemar.exact(matrix)

#save in table
mcnemar_result <- c(paste0("match_", threshold), test$estimate, test$conf.int[1], test$conf.int[2], test$p.val, length(treated_indices), length(control_indices))
odds_ratios[index,] <- mcnemar_result
index <- index + 1

#Method 2: using propensity score as an additional covariate
adj_formula <- as.formula(paste(outcome_var, " ~ " , adjustment_var, " + propensity + ", paste(covariates, collapse = " + ")))
adj_model <- glm(adj_formula, family="binomial"(link="logit"), data=edstays)


coeffs <- adj_model$coefficients
conf_ints <- data.frame(exp(confint.default(adj_model)))
conf_ints <- conf_ints@.Data
lower_confs <- conf_ints[[1]]
upper_confs <- conf_ints[[2]]
p_values <- coef(summary(adj_model))[,4]

print(summary(adj_model))

#measure treated and control here
num_treated <- length(which(edstays[[adjustment_var]]==1))
num_control <- length(which(edstays[[adjustment_var]]==0))

odds_ratios[index,] <- c("cov_adjustment", exp(coeffs[2]), lower_confs[2], upper_confs[2], p_values[2], num_treated, num_control)
index <- index + 1

#now save both sets of covariates
write.csv(odds_ratios, paste0(save_filepath,outcome_var,"_odds_ratios.csv"))
