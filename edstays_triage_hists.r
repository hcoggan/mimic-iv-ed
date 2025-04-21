library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(testit)
library(httpgd)
library(broom)
library(forcats)
library(scales)
hgd()

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(janitor)

#This file creates plots of the distribution of different patient demographic/vital characteristics
#to work out how to bin them 
#now for multiple visits


#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

edstays <- read.csv("edstays_with_binarised_complaints_and_age.csv")

print(head(edstays))
num_patients <- nrow(edstays)

edstays$admitted <- factor(edstays$admitted, levels=c(0, 1))


#for race and gender, which right now are our only time-invariant characteristics (we've dropped people changing gender as we have no way of working out if they are actually trans and there are only 383 of them)
#we need to not double-count patients

patients <- data.frame(subject_id = edstays$subject_id, race = edstays$race, gender = edstays$gender) #this is same for all patients across time
#ensure one row for each patients
patients <- patients %>% distinct(subject_id, .keep_all=TRUE)

#for race (ensure descending order)
edstays$race <- fct_infreq(edstays$race)
patients$race <- fct_infreq(patients$race)
#edstays$race <- factor(edstays$race, levels=rev(unique(edstays$race)))

#shortened_race_labels <- rev(c("American Indian", "Asian", "Black/African-American", "Hispanic/Latino", "Pacific Islander", "Other", "Unknown", "White"))
p <- ggplot(patients, aes(x = race)) + geom_bar(fill='white') + scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() +   theme_dark()
p <- p + labs(x = "Race") + labs(y="Number of patients") + labs(title = "Race of patients using the ER") + labs(caption="data from MIMIC-IV-ED, patients counted only once")
ggsave("plots/race_plot_by_patients.pdf", plot=p)

#for gender
edstays$gender <- fct_infreq(edstays$gender)
p <- ggplot(patients, aes(x = gender)) + geom_bar(width=0.5,  fill='white') + scale_fill_brewer() + theme(text = element_text(size = 20)) + coord_flip() +   theme_dark()
p <- p + labs(x = "Gender") + labs(y="Number of patients") + labs(title = "Gender of patients using the ER") + labs(caption="data from MIMIC-IV-ED, patients counted only once") 
ggsave("plots/gender_plot_by_patients.pdf", plot=p)

#everything else is time-variant, so we count ADMISSIONS, not patients

#to be able to colour by admission, we do need to start double-counting patients!

#shortened_race_labels <- rev(c("American Indian", "Asian", "Black/African-American", "Hispanic/Latino", "Pacific Islander", "Other", "Unknown", "White"))
p <- ggplot(edstays, aes(x = race, fill=admitted)) + geom_bar() + scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() +   theme_dark()
p <- p + labs(x = "Race") + labs(y="Number of admissions") + labs(title = "Patient admission by race") + labs(caption="data from MIMIC-IV-ED, patients multiple-counted")
ggsave("plots/race_plot_by_admissions.pdf", plot=p)

#for gender
edstays$gender <- fct_infreq(edstays$gender)
p <- ggplot(edstays, aes(x = gender, fill=admitted)) + geom_bar(width=0.5) + scale_fill_brewer() + theme(text = element_text(size = 20)) + coord_flip() +   theme_dark()
p <- p + labs(x = "Gender") + labs(y="Number of admissions") + labs(title = "Patient admission by gender") + labs(caption="data from MIMIC-IV-ED, patients multiple-counted") 
ggsave("plots/gender_plot_by_admissions.pdf", plot=p)


#for mode of arrival
edstays$arrival_transport <- fct_infreq(edstays$arrival_transport)
p <- ggplot(edstays, aes(x = arrival_transport, fill = admitted)) + geom_bar() + scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() +   theme_dark()
p <- p + labs(x = "Mode of arrival") + labs(y="Number of admissions") + labs(title = "Patient admission by mode of arrival") + labs(caption="data from MIMIC-IV-ED")
ggsave("plots/arrival_plot.pdf", plot=p)

#for mode of arrival for patients with NA assigned
filtered_edstays <- edstays %>% filter(acuity=="unassigned")
filtered_edstays$arrival_transport <- fct_infreq(filtered_edstays$arrival_transport)
p <- ggplot(filtered_edstays, aes(x = arrival_transport, fill = admitted)) + geom_bar() + scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() +   theme_dark()
p <- p + labs(x = "Mode of arrival") + labs(y="Number of admissions") + labs(title = "Patient admission by mode of arrival") + labs(caption="data from MIMIC-IV-ED")
ggsave("plots/arrival_plot_for_unassigned_acuity.pdf", plot=p)

#for pain-- not any obvious way to bin this
edstays$pain <- factor(edstays$pain, levels=c("none", "mild", "moderate", "severe", "UTA", "non-numeric"))
p <- ggplot(edstays, aes(x = pain, fill = admitted)) + geom_bar() + scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() +   theme_dark()
p <- p + labs(x = "Reported pain level") + labs(y="Number of admissions") + labs(title = "Patient admission by reported pain level") + labs(caption="data from MIMIC-IV-ED")
ggsave("plots/pain_plot.pdf", plot=p)

#for acuity

edstays$acuity <- factor(edstays$acuity, levels=c("low", "moderate", "high", "very_high", "unassigned"))
p <- ggplot(edstays, aes(x = acuity, fill = admitted)) + geom_bar() +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Acuity level") + labs(y="Number of admissions") + labs(title = "Patient admission by acuity score", subtitle="5=highest acuity, 1=lowest acuity, 0=no acuity assigned") + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/acuity_plot.pdf", plot=p)


#now for the fun stuff (continuous variables)
#first temperature
bin_width <- 0.2
p <- ggplot(edstays, aes(x = temperature, fill=admitted)) + geom_histogram(binwidth=bin_width) + scale_fill_brewer() + theme(text = element_text(size = 15)) + theme_dark() + coord_flip()
p <- p + labs(x = "Patient temperature") + labs(y="Number of admissions") + labs(title = "Patient admission by temperature", subtitle=paste0("to ", as.character(bin_width), " degrees")) + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/temp_plot.pdf", plot=p)

#use fever definitions given by Harvard health

edstays$fever_grade <- factor(edstays$fever_grade, levels=c("none", "low", "moderate", "high"))
p <- ggplot(edstays, aes(x = fever_grade, fill = admitted)) + geom_bar() +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Fever grade") + labs(y="Number of admissions") + labs(title = "Patient admission by fever grade", subtitle="using Harvard Health fever scale") + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/fever_plot.pdf", plot=p)


print(table(edstays$fever_grade))

#now oxygen saturation
bin_width <- 1
p <- ggplot(edstays, aes(x = o2sat, fill=admitted)) + geom_histogram(binwidth=bin_width) + scale_fill_brewer() + theme(text = element_text(size = 15)) + theme_dark() + coord_flip()
p <- p + labs(x = "Patient O2 saturation level") + labs(y="Number of admissions") + labs(title = "Patient admission by O2 saturation", subtitle=paste0("to ", as.character(bin_width), "%")) + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/o2sat_plot.pdf", plot=p)


#use Wold 2015 definitions
#edstays$o2_level <- factor(edstays$o2_level, levels=c("normal", "low", "very low"))
p <- ggplot(edstays, aes(x = o2_level, fill = admitted)) + geom_bar(width=0.5) +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "O2 saturation level") + labs(y="Number of admissions") + labs(title = "Patient admission by O2 level", subtitle="abnormal=95% or lower (Wold 2015)") + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/o2_level_plot.pdf", plot=p)

print(table(edstays$o2_level))

#now respiration rate
bin_width <- 2
p <- ggplot(edstays, aes(x = resprate, fill=admitted)) + geom_histogram(binwidth=bin_width) + scale_fill_brewer() + theme(text = element_text(size = 15)) + theme_dark() + coord_flip()
p <- p + labs(x = "Patient respiration rate") + labs(y="Number of admissions") + labs(title = "Patient admission by respiration rate", subtitle=paste0("to nearest ", as.character(bin_width), " breaths per minute")) + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/resprate_plot.pdf", plot=p)


print(table(edstays$resprate))

#use NIH definitions
edstays$resp_level <- factor(edstays$resp_level, levels=c("low", "normal", "high"))
p <- ggplot(edstays, aes(x = resp_level, fill = admitted)) + geom_bar() +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Patient respiration level") + labs(y="Number of admissions") + labs(title = "Patient admission by respiration level", subtitle="Low = less than 12 breaths per minute, high = more than 20 (NIH)") + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/resp_level_plot.pdf", plot=p)

print(table(edstays$resp_level))

#OK, now heart rate
print(table(edstays$heartrate))
bin_width <- 2
p <- ggplot(edstays, aes(x = heartrate, fill=admitted)) + geom_histogram(binwidth=bin_width) + scale_fill_brewer() + theme(text = element_text(size = 15)) + theme_dark() + coord_flip()
p <- p + labs(x = "Patient heart rate") + labs(y="Number of admissions") + labs(title = "Patient admission by heart rate", subtitle=paste0("to ", as.character(bin_width), " bpm")) + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/heart_rate_plot.pdf", plot=p)

#use NIH definitions of heart rate
edstays$heart_level <- factor(edstays$heart_level, levels=c("low", "normal", "high"))
p <- ggplot(edstays, aes(x = heart_level, fill = admitted)) + geom_bar() +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Patient heart rate level") + labs(y="Number of admissions") + labs(title = "Patient admission by heart rate level", subtitle="Low = less than 60 bpm, high = more than 100 (Nursing Skills)") + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/heart_level_plot.pdf", plot=p)

print(table(edstays$sbp))
print(table(edstays$dbp))
 
#make a density plot of blood pressure

p <- ggplot(edstays, aes(x = sbp, y=dbp)) + geom_density_2d(color="white") + scale_fill_brewer() + theme(text = element_text(size = 15)) + theme_dark() 
p <- p + labs(x = "systolic") + labs(y="diastolic") + labs(title = "Blood pressure of patients in ER")+ labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/bp_density_plot.pdf", plot=p)


#p <- ggplot(edstays_admitted, aes(x = sbp, y=dbp)) + geom_density_2d(color='White') + scale_fill_brewer() + theme(text = element_text(size = 15)) + theme_dark() 
#p <- p + labs(x = "systolic") + labs(y="diastolic") + labs(title = "Blood pressure of patients admitted")+ labs(caption="data from MIMIC-IV-ED, first visits only") 
#ggsave("plots/bp_density_plot_admitted.pdf", plot=p)

#use Johns Hopkins def of hypertension + definition of instability
edstays$bp_level <- factor(edstays$bp_level, levels=c("low", "normal", "elevated", "stage_1_high", "stage_2_high", "unstable"))
p <- ggplot(edstays, aes(x = bp_level, fill = admitted)) + geom_bar() +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Patient blood pressure level") + labs(y="Number of admissions") + labs(title = "Patient admission by blood pressure level", subtitle="Using Johns Hopkins definitions") + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/bp_level_plot.pdf", plot=p)

#POST MIMIC DATA-look at age!

bin_width <- 1
p <- ggplot(edstays, aes(x = age, fill=admitted)) + geom_histogram(binwidth=bin_width) + scale_fill_brewer() + theme(text = element_text(size = 15)) + theme_dark() + coord_flip()
p <- p + labs(x = "Patient age") + labs(y="Number of admissions") + labs(title = "Patient admission by age") + labs(caption="data from MIMIC-IV") 
ggsave("plots/age_plot.pdf", plot=p)

#bin into decades
p <- ggplot(edstays, aes(x = age_group, fill = admitted)) + geom_bar() +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Patient age") + labs(y="Number of visits") + labs(title = "Patient admission by age") + labs(caption="data from MIMIC-IV") 
ggsave("plots/age_group_plot.pdf", plot=p)


#check no pediatric patients
assert(min(edstays$age) >= 18)

#look at time of arrival
p <- ggplot(edstays, aes(x = time_of_day, fill = admitted)) + geom_bar() +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Time of arrival") + labs(y="Number of visits") + labs(title = "Patient admission by time of arrival") + labs(caption="data from MIMIC-IV") 
ggsave("plots/time_of_day_plot.pdf", plot=p)



#look at anchor-year groups
p <- ggplot(edstays, aes(x = anchor_year_group, fill = admitted)) + geom_bar() +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Year-group of admission") + labs(y="Number of visits") + labs(title = "Patient admission by year group") + labs(caption="data from MIMIC-IV") 
ggsave("plots/year_group_plot.pdf", plot=p)

adm_edstays <- read.csv("adm_edstays.csv")

#for ADMITTED PATIENTS ONLY: look at distribution of marital status
p <- ggplot(adm_edstays, aes(x = marital_status)) + geom_bar(fill='White') +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Marital status") + labs(y="Number of admissions") + labs(title = "Marital status (admitted patients only)") + labs(caption="data from MIMIC-IV") 
ggsave("plots/adm_marital_status_plot.pdf", plot=p)

#for ADMITTED PATIENTS ONLY: look at distribution of language
p <- ggplot(adm_edstays, aes(x = language)) + geom_bar(fill='White') +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Language spoken") + labs(y="Number of admissions") + labs(title = "Language spoken (admitted patients only)") + labs(caption="data from MIMIC-IV") 
ggsave("plots/adm_language_plot.pdf", plot=p)

#for ADMITTED PATIENTS ONLY: look at distribution of marital status
p <- ggplot(adm_edstays, aes(x = insurance)) + geom_bar(fill='White') +scale_fill_brewer() + theme(text = element_text(size = 15)) + coord_flip() + theme_dark()
p <- p + labs(x = "Insurance status") + labs(y="Number of admissions") + labs(title = "Insurance status (admitted patients only)") + labs(caption="data from MIMIC-IV") 
ggsave("plots/adm_insurance_status_plot.pdf", plot=p)

#for ADMITTED PATIENTS ONLY: look at distribution of wait times
bin_width <- 5
p <- ggplot(adm_edstays, aes(x = time_in_ed)) + geom_histogram(binwidth=bin_width, fill='White') + scale_fill_brewer() + theme(text = element_text(size = 15)) + theme_dark() 
p <- p + labs(x = "Wait time (minutes)") + labs(y="Number of admissions") + labs(title = "Frequency of wait times for admitted patients", subtitle=paste0("to ", as.character(bin_width), " minutes")) + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/adm_wait_time_plot.pdf", plot=p)

print("Distribution of wait times for admitted patients")
print(quantile(adm_edstays$time_in_ed, probs=c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9)))


#for DISCHARGED PATIENTS ONLY: look at distribution of wait times
edstays <- edstays %>% filter(disposition=="HOME")
p <- ggplot(adm_edstays, aes(x = time_in_ed)) + geom_histogram(binwidth=bin_width, fill="White") + scale_fill_brewer() + theme(text = element_text(size = 15)) + theme_dark()
p <- p + labs(x = "Wait time (minutes)") + labs(y="Number of visits") + labs(title = "Frequency of wait times for discharged patients", subtitle=paste0("to ", as.character(bin_width), " minutes")) + labs(caption="data from MIMIC-IV-ED") 
ggsave("plots/discharged_wait_time_plot.pdf", plot=p)

print("Distribution of wait times for discharged patients")
print(quantile(edstays$time_in_ed, probs=c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9)))