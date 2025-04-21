library(lubridate)
library(ggsurvfit)
#library(gtsummary)
library(tidycmprsk)
library(testit)
library(httpgd)
library(purrr)
hgd()

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(rlang)

#This file converts times in edstays to 'minutes since baseline date' 
#adds new columns to edstays indicating how many previous visits there have been within certain timeframes
#it also records the time of day at which patients are admitted
#then it saves

#setwd("/Users/helenacoggan/Documents/MIMIC-IV-ED/")

edstays <- read.csv("mimic-iv-ed-2.2/ed/edstays.csv")


#any NA values before conversion?
print("before conversion")
print(which(is.na(edstays$intime)))
print(which(is.na(edstays$outtime))) #all are machine readable before conversion


#convert dates to a format R can manipulate
outtimes <- as.POSIXct(edstays$outtime, format = "%Y-%m-%d %H:%M:%S")
intimes <- as.POSIXct(edstays$intime, format = "%Y-%m-%d %H:%M:%S")

#record time of arrival
#add column to dataframe
zero_vector <- c(rep(0, nrow(edstays)))
edstays$time_of_day <- zero_vector

#now classify each arrival time
hour_of_arrival <- hour(edstays$intime)
indices <- which((hour_of_arrival>=0) & (hour_of_arrival<6))
edstays$time_of_day[indices] <- 'small hours'
indices <- which((hour_of_arrival>=6) & (hour_of_arrival<12))
edstays$time_of_day[indices] <- 'morning'
indices <- which((hour_of_arrival>=12) & (hour_of_arrival<18))
edstays$time_of_day[indices] <- 'afternoon'
indices <- which((hour_of_arrival>=18) & (hour_of_arrival<24))
edstays$time_of_day[indices] <- 'evening'


#record year of arrival
edstays$year_of_arrival <- year(edstays$intime)

#set a baseline date
baseline_date <- as.POSIXct("2009-12-31 00:00:00", format = "%Y-%m-%d %H:%M:%S") #24 hours before the start of the 2010s, since visits start in 2011

#for purposes of comparison with date of death, keep intime as a column and rename it
edstays$intime_as_date <- edstays$intime


edstays$outtime <- difftime(outtimes, baseline_date, units="mins")
edstays$intime <- difftime(intimes, baseline_date, units="mins") #measure all times in minutes since some ref point


#any NA values after conversion?
print("after conversion")
print(length(which(is.na(edstays$intime))))
print(length(which(is.na(edstays$outtime)))) #a small number (26 and 45 respectively) aren't legible to the converter; we delete them

edstays <- edstays %>%
  filter(!is.na(intime) & !is.na(outtime), outtime > intime)


#


# edstays <- edstays %>% 
#     group_by(subject_id) %>%
#     mutate(previous_visits = row_number(intime) - 1) %>% #0 if 1st visit, etc
#     ungroup()

# #now look at admissions
# adm <- edstays %>% 
#     filter(disposition=="ADMITTED") %>%
#     group_by(subject_id) %>%
#     mutate(previous_admissions = row_number(intime) - 1) %>% #0 if 1st admission, etc
#     ungroup()

# #now join this up
# adm <- subset(adm, select=c("stay_id", "previous_admissions"))
# edstays <- edstays %>% left_join(adm, by="stay_id")

# edstays$previous_admissions[which(is.na(edstays$previous_admissions))] <- 0 #if no value, then no previous visits

#save modified edstays 

#add flag for number of visits per patient, and filter patients with more than 1 visit
edstays <- edstays %>% add_count(subject_id)
df <- edstays %>% filter(n > 1)

window <- 30*24*60 #1 month in mins


#create function
prev_visit <- function(i){
    id <- df$subject_id[i]
    t <- df$intime[i]
    time_threshold <- t - window
    indices <- which(df$subject_id==id & df$intime >= time_threshold & df$intime < t)
    prev_visits <- df[indices,]
    if (length(indices)> 0) {
      if (any(prev_visits$disposition=="ADMITTED")) {
        return("recent_admission") #there has been a recent admission
      } else {
        return("recent_visit_without_admission") #there has been a recent visit without admission
      }
    } else {
      return("no_recent_visits")
    }
}

df$recency <- purrr::map_vec(seq(1:nrow(df)), prev_visit, .progress=TRUE)
#print(recency)

write.csv(df, "edstays_mult_visit_time_window_recorded.csv")


#df <- read.csv("edstays_mult_visit_time_window_recorded.csv")

#limit to key columns
df <- df %>% select(c("stay_id", "recency"))

#now join back to edstays
edstays <- edstays %>% left_join(df, by='stay_id')
edstays$recency[is.na(edstays$recency)] <- "no_recent_visits" #for patients with no recent visits
cols <- colnames(edstays)
cols <- cols[!(cols == "n")]
edstays <- edstays[cols]

edstays <- as.data.frame(edstays)

write.csv(edstays, "edstays_mult_visit_with_timeframes_without_vitals.csv")


