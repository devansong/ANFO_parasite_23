##################################################################
##################################################################
######ANFO parasite data analysis for Memphis Zoo ################
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

source("~MRFcov_ID.R")
source("~cv_MRF_diag_rep_ID.R")
source("~bootstrap_MRF_ID.R")

source("/Users/annedevan-song_1/Dropbox/Publications_Work/ANFO_Parasites/Code/MRFcov_ID.R")
source("/Users/annedevan-song_1/Dropbox/Publications_Work/ANFO_Parasites/Code/bootstrap_MRF_ID.R")
source("/Users/annedevan-song_1/Dropbox/Publications_Work/ANFO_Parasites/Code/cv_MRF_diag_rep_ID.R")

setwd("~/Dropbox/Publications_Work/ANFO_Parasites/Data") 



df <- read.csv("For_CRF_Poisson.csv")
##This is  to remove outliers for a 2nd round  df <- subset(df, Animal.ID != "ANFO_193" & Animal.ID!= "ANFO_29")
df <- df[, -6]

for(i in c(6:ncol(df))){
  df[,i] <- scale(as.numeric(df [,i]))
}

str(df)
#Make sure all parasites are integers 
#Correct it if it is not.

analysis.data <- df[, -1]


for(i in c(5:ncol(analysis.data))){
  analysis.data[,i] <- scale(as.numeric(analysis.data[,i]))
}

str(df)
#Make sure all parasites columns are sstructured as integer 
#Correct it if it is not 

MRF_fit <- MRFcov_ID(data = analysis.data, n_nodes = 4, n_cores = 4,
                     family = 'poisson', id_data = df, bootstrap = FALSE)
MRF_fit$key_coefs

evaluate <- cv_MRF_diag_rep_ID(data = analysis.data, n_nodes = 4,
                               n_cores = 4, family = 'poisson',id_data = df, plot = F, compare_null = T)
evaluate_nocov <- cv_MRF_diag_rep_ID(data = analysis.data[1:4], n_nodes = 4,
                                     n_cores = 4, family = 'poisson',id_data = df, plot = F, compare_null = T)

quantile(evaluate$mean_sensitivity[evaluate$model == 'Spatial MRF'], probs = c(0.05, 0.95)) 
mean(evaluate$mean_sensitivity[evaluate$model == 'Spatial MRF']) 
quantile(evaluate$mean_specificity[evaluate$model == 'Spatial MRF'], probs = c(0.05, 0.95))
mean(evaluate$mean_specificity[evaluate$model == 'Spatial MRF']) 
quantile(evaluate$mean_tot_pred[evaluate$model == 'Spatial MRF'], probs = c(0.05, 0.95))
mean(evaluate$mean_tot_pred[evaluate$model == 'Spatial MRF']) 
###
quantile(evaluate_nocov$mean_sensitivity[evaluate_nocov$model == 'Spatial MRF'], probs = c(0.05, 0.95)) 
mean(evaluate_nocov$mean_sensitivity[evaluate_nocov$model == 'Spatial MRF']) 
quantile(evaluate_nocov$mean_specificity[evaluate_nocov$model == 'Spatial MRF'], probs = c(0.05, 0.95)) 
mean(evaluate_nocov$mean_specificity[evaluate_nocov$model == 'Spatial MRF'])
quantile(evaluate_nocov$mean_tot_pred[evaluate_nocov$model == 'Spatial MRF'], probs = c(0.05, 0.95)) 
mean(evaluate_nocov$mean_tot_pred[evaluate_nocov$model == 'Spatial MRF'])


#####
booted_CRF_all <- bootstrap_MRF_ID(data = analysis.data, n_bootstraps = 100,
                                   n_nodes = 4, n_cores = 1, family = 'poisson',
                                   id_data = df, sample_seed = 122696)

#Parameters that should not be changed for this line are n_nodes 
#####
booted_CRF_all$mean_key_coefs #This should closely if not completely match table 2 in manuscript. 
newdf <- as.data.frame(booted_CRF_all$mean_key_coefs) 






