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
library(MRFcov)
library(parallel)
library(igraph)
library(PerformanceAnalytics)
library(fastDummies)
library(ggpubr)
library(gridExtra)
library(cowplot)

source("~MRFcov_ID.R") #This R document is uploaded to Github
source("~cv_MRF_diag_rep_ID.R") #This R document is uploaded to Github
source("~bootstrap_MRF_ID.R") #This R document is uploaded to Github


df <- read.csv("For_CRF_Binomial.csv")
df <- df[, -6]

for(i in c(6:ncol(df))){
  df[,i] <- scale(as.numeric(df [,i]))
}

str(df)
#Make sure all parasites are binary (0 and 1)
#Correct it if it is not. Node data should be binary 

analysis.data <- df[, -1]

####################################################################################################
MRF_fit <- MRFcov_ID(data = analysis.data, n_nodes = 4, n_cores = 4,
                     family = 'binomial', id_data = df, bootstrap = FALSE)

plotMRF_hm(MRF_mod = MRF_fit, main = 'MRF (no covariates)_POISSON', 
           node_names = c("Strongyle","Coccidia","Tapeworm","Ascarid"))

MRF_fit$key_coefs

evaluate <- cv_MRF_diag_rep_ID(data = analysis.data, n_nodes = 4,
                               n_cores = 4, family = 'binomial',id_data = df, plot = F, compare_null = T)
evaluate_nocov <- cv_MRF_diag_rep_ID(data = analysis.data[1:4], n_nodes = 4,
                                     n_cores = 4, family = 'binomial',id_data = df, plot = F, compare_null = T)
#####
quantile(evaluate$mean_sensitivity[evaluate$model == 'Spatial MRF'], probs = c(0.05, 0.95)) 
mean(evaluate$mean_sensitivity[evaluate$model == 'Spatial MRF']) 
quantile(evaluate$mean_specificity[evaluate$model == 'Spatial MRF'], probs = c(0.05, 0.95)) 
mean(evaluate$mean_specificity[evaluate$model == 'Spatial MRF']) 
quantile(evaluate$mean_tot_pred[evaluate$model == 'Spatial MRF'], probs = c(0.05, 0.95)) 
mean(evaluate$mean_tot_pred[evaluate$model == 'Spatial MRF']) 
######
quantile(evaluate_nocov$mean_sensitivity[evaluate_nocov$model == 'Spatial MRF'], probs = c(0.05, 0.95))  
mean(evaluate_nocov$mean_sensitivity[evaluate_nocov$model == 'Spatial MRF']) 
quantile(evaluate_nocov$mean_specificity[evaluate_nocov$model == 'Spatial MRF'], probs = c(0.05, 0.95))
mean(evaluate_nocov$mean_specificity[evaluate_nocov$model == 'Spatial MRF']) 
quantile(evaluate_nocov$mean_tot_pred[evaluate_nocov$model == 'Spatial MRF'], probs = c(0.05, 0.95)) 
mean(evaluate_nocov$mean_tot_pred[evaluate_nocov$model == 'Spatial MRF']) 
#####
booted_CRF_all <- bootstrap_MRF_ID(data = analysis.data, n_bootstraps = 100,
                                   n_nodes = 4, n_cores = 1, family = 'binomial',
                                   id_data = df, sample_seed = 122696)

#Experiment with changing number of cores if this takes too long 
#Parameters that should not be changed for this line are n_nodes 

#######
booted_CRF_all$mean_key_coefs #This should closely if not completely match table 1 in manuscript. 
newdf <- as.data.frame(booted_CRF_all$mean_key_coefs) 

