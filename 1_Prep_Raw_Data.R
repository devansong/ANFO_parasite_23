##################################################################
##################################################################
######ANFO parasite data analysis for Memphis Zoo project#########
####Code by Anne Devan-Song. Bend, OR. 2023#######################
##################################################################
##################################################################


rm(list=ls())
graphics.off()

library(dplyr)
library(tidyverse)
library(dplyr)
library(MRFcov)
library(parallel)
library(igraph)
library(PerformanceAnalytics)
library(fastDummies)
library(ggpubr)
library(gridExtra)
library(cowplot)

dat <- read.csv("data_MERGED.csv")

df <- subset(dat, SVL>0)
temp <- as.data.frame(table(df$Animal.ID)) #Calculate recaps
mean(dat$SVL, na.rm=TRUE)
(sd(dat$SVL, na.rm=TRUE))/(sqrt(97))
mean(dat$Mass, na.rm=TRUE)
(sd(dat$Mass, na.rm=TRUE))/(sqrt(97))

table(dat$Site)
dat <- dat[complete.cases(dat),]
df <- dat
subdf <- dummy_cols(df, select_columns = 'Site') #create dummy variables for population 
df <-subdf[, !(colnames(subdf) %in% c("Site", "Site_Village Creek"))] 
subdf <- dummy_cols(df, select_columns = 'Sex') #create dummy variables for population 
subdf <- subset(subdf, Sex !="unk")
df <-subdf[, !(colnames(subdf) %in% c("Sex", "Sex_unk"))] #remove site column

mite <- as.data.frame(df$Mite.LPG)
mite[mite > 0] <- 1 
colnames(mite) <- "Mite.PA"
df$Mite.PA <- mite$Mite.PA
df <- df[, -c(3, 11)]

for(i in c(2:5)){
  df[,i] <- as.integer(df [,i]) #make binary presence/absence as numeric 
}
write.csv(df, file="For_CRF_Poisson.csv", row.names=FALSE)

paronly <- df[, c("Strongyle.EPG",
                  "Coccidia.EPG", 
                  "Cestode.EPG", 
                  "Ascarid.EPG")]
paronly[paronly > 0] <- 1 
newdf <- cbind(df[, 1], paronly, df[, 6:16])

names(newdf)[names(newdf) == "df[, 1]"] <- 'Animal.ID'
write.csv(newdf, file="For_CRF_Binomial.csv", row.names=FALSE)


