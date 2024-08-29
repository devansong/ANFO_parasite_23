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

source("~MRFcov_ID.R")
source("~cv_MRF_diag_rep_ID.R")
source("~bootstrap_MRF_ID.R")

df <- read.csv("For_CRF_Poisson.csv")
analysis.data <- df[, c("Strongyle.EPG",
                        "Coccidia.EPG",
                        "Cestode.EPG", 
                        "Ascarid.EPG", 
                        "Fecal.Mass.g.",  
                        "SVL",             
                        "BCI",            
                        "Site_Park",   
                        "Site_Forest", 
                        "Site_Zoo",        
                        "Sex_F",           
                        "Sex_M", 
                        "Mite_PA")]

for(i in c(5:ncol(analysis.data))){
  analysis.data[,i] <- scale(as.numeric(analysis.data[,i]))
}

str(df)
#Make sure all parasites columns are sstructured as integer 
#Correct it if it is not 

MRF_fit <- MRFcov_ID(data = analysis.data, n_nodes = 4, n_cores = 4,
                     family = 'poisson', id_data = df, bootstrap = FALSE)
MRF_fit$key_coefs
MRF_fit$key_coefs$Strongyle
MRF_fit$key_coefs$Coccidia
MRF_fit$key_coefs$Tapeworm
MRF_fit$key_coefs$Ascarid

evaluate <- cv_MRF_diag_rep_ID(data = analysis.data, n_nodes = 4,
                               n_cores = 4, family = 'poisson',id_data = df, plot = F, compare_null = T)
evaluate_nocov <- cv_MRF_diag_rep_ID(data = analysis.data[1:5], n_nodes = 4,
                                     n_cores = 4, family = 'binomial',id_data = df, plot = F, compare_null = T)
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


#Experiment with changing number of cores if this takes too long 
#Parameters that should not be changed for this line are n_nodes 
#####

str <- as.data.frame(booted_CRF_all$mean_key_coefs$Strongyle_egg)
str$parasite <- "Strongyle"
coc <- as.data.frame(booted_CRF_all$mean_key_coefs$Coccidia_egg)
coc$parasite <- "Coccidia"
tap <- as.data.frame(booted_CRF_all$mean_key_coefs$Tapeworm_egg)
tap$parasite <- "Tapeworm"
asc <- as.data.frame(booted_CRF_all$mean_key_coefs$Ascarid_egg)
asc$parasite <- "Ascarid"
summarydf <- rbind(str,coc, tap, asc)