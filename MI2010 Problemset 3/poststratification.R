require("ipumsr")
library(haven)
library(tidyverse)
library(janitor)

# Read in the raw data.

setwd("C:/Users/sa/Documents/problemset3")
# Read in the raw data 
data <- read.table("usa_00001.dat")
# Add the labels

data <- labelled::to_factor(data)



head(data)

#Load data

ddi <- read_ipums_ddi("usa_00001.xml")
ddi <- read_ipums_micro(ddi)

#Preview
head(ddi)

#Clean

ddi$age <- as.integer(ddi$AGE)
ddi$race <- as.integer(ddi$RACE)
ddi$edu <- as.integer(ddi$EDUC)
ddi$income <- as.integer(ddi$FTOTINC)

write_csv(ddi, "census_data.csv")



         