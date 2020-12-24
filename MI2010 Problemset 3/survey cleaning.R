library(haven)
library(tidyverse)
library(ipumsr)
library(janitor)

setwd("C:/Users/sa/Documents/problemset3")
# Read in the raw data 
raw_data <- read_dta("ns20200625/ns20200625.dta")
# Add the labels
raw_data <- labelled::to_factor(raw_data)
# Just keep some variables
reduced_data <- 
  raw_data %>% 
  select(interest,
         registration,
         vote_2016,
         vote_intention,
         vote_2020,
         ideo5,
         employment,
         foreign_born,
         gender,
         census_region,
         hispanic,
         race_ethnicity,
         household_income,
         education,
         state,
         congress_district,
         age)

# remove voters who arent registered from the dataset

reduced_data <- reduced_data[reduced_data$registration == "Registered" 
                             & reduced_data$vote_intention != "No, I will not vote but I am eligible", ]



#### What else???? ####
# Maybe make some age-groups?
# Maybe check the values?
# Is vote a binary? If not, what are you going to do?

reduced_data<-
  reduced_data %>%
  mutate(vote_trump = 
           ifelse(vote_2020=="Donald Trump", 1, 0))

reduced_data<-
  reduced_data %>%
  mutate(race = 
           ifelse(race_ethnicity=="White", 1, 0))

reduced_data<-
  reduced_data %>%
  rowwise() %>% 
  mutate(ideal = str_remove(string = ideo5, pattern = "\\ .*")) %>% 
  mutate(ideal = case_when(
    ideo5=="Very Liberal" ~ 1,
    ideo5=="Liberal" ~ 2,
    ideo5=="Moderate" ~ 3,
    ideo5=="Conservative" ~ 4,
    ideo5=="Very Conservative" ~ 5,
    ideo5=="Not Sure" ~ 0
  )) 

reduced_data<-
  reduced_data %>%
  rowwise() %>% 
  mutate(edu = str_remove(string = education, pattern = "\\ .*")) %>% 
  mutate(edu = case_when(
    education == "3rd Grade or less" ~ 1,
    education == "Middle School - Grades 4 - 8" ~ 2,
    education == "Completed some high school" ~ 3,
    education == "High school graduate" ~ 4,
    education == "Other post high school vocational training" ~ 5,
    education == "Completed some college, but no degree" ~ 6,
    education == "Associate Degree" ~ 7,
    education == "College Degree (such as B.A., B.S.)" ~ 8,
    education == "Completed some graduate, but no degree" ~ 9,
    education == "Masters degree" ~ 10,
    education == "Doctorate degree" ~ 11,
    education == "NA" ~ 0
  )) 

reduced_data<-
  reduced_data %>%
  rowwise() %>% 
  mutate(income = str_remove(string = household_income, pattern = "\\ .*")) %>% 
  mutate(income = case_when(
    household_income=="$15,000 to $19,999" ~ 15000,
    household_income=="$20,000 to $24,999" ~ 20000,
    household_income=="$25,000 to $29,999" ~ 25000,
    household_income=="$30,000 to $34,999" ~ 30000,
    household_income=="$35,000 to $39,999" ~ 35000,
    household_income=="$40,000 to $44,999" ~ 40000,
    household_income=="$45,000 to $49,999" ~ 45000,
    household_income=="$50,000 to $54,999" ~ 50000,
    household_income=="$55,000 to $59,999" ~ 55000,
    household_income=="$60,000 to $64,999" ~ 60000,
    household_income=="$65,000 to $69,999" ~ 65000,
    household_income=="$70,000 to $74,999" ~ 70000,
    household_income=="$75,000 to $79,999" ~ 75000,
    household_income=="$80,000 to $84,999" ~ 80000,
    household_income=="$85,000 to $89,999" ~ 85000,
    household_income=="$90,000 to $94,999" ~ 90000,
    household_income=="$95,000 to $99,999" ~ 95000,
    household_income=="$100,000 to $124,999" ~ 100000,
    household_income=="$125,000 to $149,999" ~ 125000,
    household_income=="$150,000 to $199,999" ~ 150000,
    household_income=="$200,000 to $249,999" ~ 200000,
    household_income=="$250,000 and above" ~ 100000,
    household_income=="NA" ~ 0
  )) 


#take a look at the datatypes, etc
head(reduced_data)

# Saving the survey/sample data as a csv file in my
# working directory
write_csv(reduced_data, "survey_data.csv")

logmodel <- glm( vote_trump ~ income + age + edu + ideal + race, family = binomial, 
             model = TRUE, method = "glm.fit", data = reduced_data)

summary(logmodel)

