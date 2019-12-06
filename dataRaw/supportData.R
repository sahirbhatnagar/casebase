######################################
# R Source code file for creating dataset to be included in the sail package
# data taken from:
# http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/support.tsv
# which is from the Study to Understand Prognoses Preferences Outcomes and Risks of Treatment (SUPPORT)
# http://biostat.mc.vanderbilt.edu/wiki/Main/SupportDesc
# Code modified from https://github.com/sahirbhatnagar/sail/blob/master/data-raw/OASIS_data.R
# Author: Jesse Islam
# Created: December 6, 2019
# Updated:
#####################################

## Read in the SUPPORT data

if (!require("pacman")) install.packages("pacman")
pacman::p_load(magrittr)
pacman::p_load(dplyr)
pacman::p_load(usethis)
pacman::p_load(readr)
pacman::p_load(mice)
support <- read_tsv("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/support.tsv")
#we will impute the data:


#impute without death, to avoid correlations
sup.impute <-mice(support[,-c(2)], m = 5, maxit = 10,printFlag=FALSE)
#take a complete dataset, Only one is being made and used here.
#The logic being we only want to compare between models.
#Any Bias introduced is introduced to all models.
completeData<-complete(sup.impute,1)
#put death back into the complete dataset, as it was removed.
completeData$death<-support$death
#make everything into factors and remove observations without
#imputed values.
completeData<-na.omit(as.data.frame(unclass(completeData)))
support <- completeData

usethis::use_data(support, overwrite = TRUE)
