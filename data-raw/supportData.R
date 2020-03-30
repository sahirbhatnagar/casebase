######################################
# R Source code file for creating dataset to be included in the sail package
# data taken from:
# http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/support2csv.zip
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

#DATA retrieved from http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/support2csv.zip
ew <- read_csv("data-raw/support2.csv")
data<-read.csv('data-raw/support2.csv', col.names=colnames(ew), fill=TRUE, header=TRUE)
rm(ew)
#manual imputation
data$alb[is.na(data$alb)]=3.5
data$pafi[is.na(data$pafi)]=333.3
data$bili[is.na(data$bili)]=1.01
data$crea[is.na(data$crea)]=1.01
data$bun[is.na(data$bun)]=6.51
data$wblc[is.na(data$wblc)]=9
data$urine[is.na(data$urine)]=2502

# automated imputation
previousModelPositions<-c(18,19,20,21,26,27,28,29)

supportMain<-data[,-previousModelPositions]
supportPrevious<-data[,previousModelPositions]
#keep 1,2 after imputation from previous usable models.

#missingPattern<-md.pattern(supportMain)
supportMainImpute <-mice(supportMain[,-c(4,11,38,13,14,15)], m = 1,printFlag=FALSE) #impute while removing ordinal variables and other response variables.
imputedData<-cbind(complete(supportMainImpute,1),supportPrevious[,c(1,2)])
#adls is missing 33%, and is removed accordingly and hospitaldeath
#md.pattern(completeData)
completeData<-na.omit(imputedData[,-c(32)])
support <- completeData

usethis::use_data(support, overwrite = TRUE)
