bootstrap_MRF_ID <- function(data, n_bootstraps, sample_seed, symmetrise,
                          n_nodes, n_cores, n_covariates, family,
                          sample_prop,
                          spatial = TRUE,
                          id_data){
  
  #### Specify default parameter values and initiate warnings #####
  if(!(family %in% c('gaussian', 'poisson', 'binomial')))
    stop('Please select one of the three family options:
         "gaussian", "poisson", "binomial"')
  
  if(missing(symmetrise)){
    symmetrise <- 'mean'
  }
  
  if(!(symmetrise %in% c('min', 'max', 'mean')))
    stop('Please select one of the three options for symmetrising coefficients:
         "min", "max", "mean"')
  
  if(missing(n_cores)) {
    n_cores <- 1
  } else {
    if(sign(n_cores) != 1){
      stop('Please provide a positive integer for n_cores')
    } else{
      if(sfsmisc::is.whole(n_cores) == FALSE){
        stop('Please provide a positive integer for n_cores')
      }
    }
  }
  
  if(missing(sample_seed)){
    set.seed(ceiling(runif(1, 0, 100000)))
  } else {
    set.seed(sample_seed)
  }
  
  if(missing(n_bootstraps)) {
    n_bootstraps <- 10
  } else {
    if(sign(n_bootstraps) == 1){
      #Make sure n_bootstraps is a positive integer
      n_bootstraps = ceiling(n_bootstraps)
    } else {
      stop('Please provide a positive integer for n_bootstraps')
    }
  }
  
  
  lambda1 <- 1
  
  #### Basic checks on data arguments ####
  if(missing(n_nodes)) {
    warning('n_nodes not specified. using ncol(data) as default, assuming no covariates',
            call. = FALSE)
    n_nodes <- ncol(data)
    n_covariates <- 0
  } else {
    if(sign(n_nodes) != 1){
      stop('Please provide a positive integer for n_nodes')
    } else {
      if(sfsmisc::is.whole(n_nodes) == FALSE){
        stop('Please provide a positive integer for n_nodes')
      }
    }
  }
  
  if(missing(n_covariates)){
    n_covariates <- ncol(data) - n_nodes
  } else {
    if(sign(n_covariates) != 1){
      stop('Please provide a positive integer for n_covariates')
    } else {
      if(sfsmisc::is.whole(n_covariates) == FALSE){
        stop('Please provide a positive integer for n_covariates')
      }
    }
  }
  
  if(n_covariates > 0){
    if(any(is.na(data[,(n_nodes + 1):ncol(data)]))){
      warning('NAs detected in covariate columns. These will be imputed from rnorm(mean=0,sd=1)',
              call. = FALSE)
      nas_present <- TRUE
    } else {
      nas_present <- FALSE
    }
  } else {
    nas_present <- FALSE
  }
  
  if(nrow(data) < 2){
    stop('The data must have at least 2 rows')
  }
  
  if(any(!(apply(data, 2, class) %in% c('numeric', 'integer')))){
    stop('Only integer and numeric values permitted')
  }
  
  if(is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  #### Use paranormal transformation for Poisson variables ####
  if(family == 'poisson'){
    cat('Poisson variables will be transformed using a nonparanormal...\n')
    
    #square_root_mean = function(x) {sqrt(mean(x ^ 2))}
    #poiss_sc_factors <- apply(data[, 1:n_nodes], 2, square_root_mean)
    #data[, 1:n_nodes] <- apply(data[, 1:n_nodes], 2,
    #                           function(x) x / square_root_mean(x))
    
    # Function to estimate parameters of a nb distribution
    nb_params = function(x){
      MASS::fitdistr(x, densfun = "negative binomial")$estimate
    }
    
    # Function to estimate parameters of a poisson distribution
    poiss_params = function(x){
      MASS::fitdistr(x, densfun = "poisson")$estimate
    }
    
    # Function to transform counts using nonparanormal
    paranorm = function(x){
      ranks <- rank(log2(x + 0.01))
      stats::qnorm(ranks / (length(x) + 1))
    }
    
    # Calculate raw parameters
    suppressWarnings(poiss_sc_factors <- try(apply(data[, 1:n_nodes],
                                                   2, nb_params), silent = TRUE))
    
    if(inherits(poiss_sc_factors, 'try-error')){
      suppressWarnings(poiss_sc_factors <- apply(data[, 1:n_nodes],
                                                 2, poiss_params))
    }
    
    data[, 1:n_nodes] <- apply(data[, 1:n_nodes], 2, paranorm)
    family <- 'gaussian'
    return_poisson <- TRUE
    
  } else {
    return_poisson <- FALSE
  }
  
  #### Function to randomly sample rows for each bootstrap replicate ####
  if(missing(sample_prop)){
    sample_prop <- 1
  }
  
  if(sample_prop < 0.1 || sample_prop > 1){
    stop('sample_prop must be a proportion ranging from 0.1 to 1')
  }
  
  shuffle_rows <- function(x){
    dplyr::sample_n(x, nrow(x) * sample_prop, FALSE)
  }
  
  #### Function to impute covariate NAs from normal distribution (mean = 0; sd = 1) ####
  impute_nas <- function(empty){
    data[is.na(data)] <- sample(rnorm(sum(is.na(data)),
                                      mean = 0, sd = 1),
                                replace = FALSE)
    data <- data
  }
  
  #### Function to count proportions of non-zero coefficients ####
  countzero <- function(data, x, y){
    bs.unlist <- data %>% purrr::map('direct_coefs')
    estimatesinxy <- unlist(lapply(bs.unlist, '[', x, y))
    zeron <- length(which(estimatesinxy == 0))
    
    #correct for finite sampling
    ((.0001 * length(data)) + zeron) / ((.0001 * length(data)) + length(data))
  }
  
  #### Run MRFcov across iterations ####
  #Needs to be a sequence for running models repeatedly. Uses n_cores to define the length
  lambda1_seq <- seq_len(n_bootstraps / n_cores)
  
  
  #### Specify number of bootstrap replicates for each processing core to run ####
  #Set the total number of repititions
  total_reps <- n_bootstraps
  n_bootstraps <- suppressWarnings(lapply(split(seq_len(n_bootstraps),
                                                lambda1_seq),
                                          length))
  
  #### If n_cores > 1, check parallel library loading ####
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)
    
    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)
    
    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {
      
      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/purrr$", "", system.file(package = "purrr"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))
      
      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)
      
      if(class(test_load2) == "try-error"){
        
        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(penalized)), silent = TRUE)
        
        if(class(test_load3) == "try-error"){
          
          #Give up and use lapply instead!
          parallel_compliant <- FALSE
          stopCluster(cl)
          
        } else {
          parallel_compliant <- TRUE
        }
        
      } else {
        parallel_compliant <- TRUE
      }
      
    } else {
      parallel_compliant <- TRUE
    }
  } else {
    #If n_cores = 1, set parallel_compliant to FALSE
    parallel_compliant <- FALSE
  }
  
  #### If parallel support confirmed and n_cores > 1, proceed with parLapply ####
  if(parallel_compliant){
    cat('Fitting bootstrap_MRF models in parallel using', n_cores, 'cores... \n')
    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('lambda1_seq',
                          'symmetrise', 'n_nodes', 'n_bootstraps',
                          'n_covariates', 'family', 'nas_present',
                          'spatial', 'data','id_data'),
                  envir = environment())
    
    if(spatial){
      clusterExport(NULL, c('MRFcov_ID'), envir = environment())
      
    }
    
    #Export necessary functions to each cluster
    clusterExport(NULL, c('MRFcov', 'prep_MRF_covariates'))
    clusterExport(NULL, 'countzero', envir = environment())
    clusterExport(NULL, 'shuffle_rows', envir = environment())
    clusterExport(NULL, 'sample_prop', envir = environment())
    
    #Export necessary libraries
    clusterEvalQ(cl, library(purrr))
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(data.table))
    
    lambda_results <- pbapply::pblapply(lambda1_seq, function(l) {
      booted_mrfs <- lapply(seq_len(n_bootstraps[[l]]), function(x) {
        sample_data <- sample(seq_len(100), 1)
        
        if(spatial){
          ID <- as.data.frame(id_data$Animal.ID)
          colnames(ID) <- 'Animal.ID'
          booted_data <- shuffle_rows(cbind(data, ID))
          sample_data <- booted_data[, 1:(ncol(booted_data) - ncol(ID))]
          
          if(nas_present){
            sample_data <- impute_nas(sample_data)
          }
          
          sample_data <- prep_MRF_covariates(sample_data, n_nodes)
          sample_id <- booted_data[, ((ncol(booted_data) + 1) -
                                            ncol(ID)):ncol(booted_data)]
          sample_id <- as.data.frame(sample_id)
          colnames(sample_id) <- 'Animal.ID'
          sample_id$Animal.ID <- as.factor(sample_id$Animal.ID)
          
          mod <- suppressWarnings(MRFcov_ID(data = sample_data,
                                                 symmetrise = symmetrise,
                                                 id_data = sample_id,
                                                 n_nodes = n_nodes,
                                                 n_cores = 1,
                                                 prep_covariates = FALSE,
                                                 n_covariates = n_covariates,
                                                 family = family,
                                                 bootstrap = TRUE))
          
        } else {
          sample_data <- shuffle_rows(data)
          
          if(nas_present){
            sample_data <- impute_nas(booted_data)
          }
          
          sample_data <- prep_MRF_covariates(sample_data, n_nodes)
          mod <- suppressWarnings(MRFcov(data = sample_data,
                                         symmetrise = symmetrise,
                                         n_nodes = n_nodes,
                                         n_cores = 1,
                                         prep_covariates = FALSE,
                                         n_covariates = n_covariates,
                                         family = family,
                                         bootstrap = TRUE))
        }
        
        list(direct_coefs = as.matrix(mod$direct_coefs),
             indirect_coefs = as.matrix(mod$indirect_coefs))
      })
      
      #Gather direct effect estimates from all bootstrap samples
      direct_coef_list <- booted_mrfs %>%
        purrr::map('direct_coefs')
      
      #Calculate mean coefficient estimates across bootstrap samples
      direct_coef_means <- apply(array(unlist(direct_coef_list),
                                       c(nrow(booted_mrfs[[1]]$direct_coefs),
                                         ncol(booted_mrfs[[1]]$direct_coefs),
                                         length(booted_mrfs))), c(1, 2), mean)
      
      rownames(direct_coef_means) <- rownames(booted_mrfs[[1]]$direct_coefs)
      colnames(direct_coef_means) <- colnames(booted_mrfs[[1]]$direct_coefs)
      
      #Calculate proportion of bootstrap models in which each cofficient is non-zero
      n_total_covariates <- ncol(booted_mrfs[[1]]$direct_coefs)
      
      prop_covs_retained <- matrix(0, n_nodes, n_total_covariates)
      for(i in seq_len(n_nodes)){
        for(j in seq_len(n_total_covariates)){
          prop_covs_retained[i, j] <- 1 - countzero(booted_mrfs, i, j)
        }
      }
      rownames(prop_covs_retained) <- rownames(booted_mrfs[[1]]$direct_coefs)
      colnames(prop_covs_retained) <- colnames(booted_mrfs[[1]]$direct_coefs)
      
      # Remove spatial splines from the key covariates list
      if(spatial){
        prop_covs_retained <- prop_covs_retained[, !grepl('Spatial', colnames(prop_covs_retained))]
        direct_coef_means_nonspat <- direct_coef_means[, !grepl('Spatial', colnames(direct_coef_means))]
      } else {
        direct_coef_means_nonspat <- direct_coef_means
      }
      
      #Create list of covariates that are not zero for each node in data
      key_covariates <- lapply(seq_len(n_nodes), function(x){
        cov_prop_retained <- apply(data.frame(prop_covs_retained)[x,], 2,
                                   function(j) ifelse(j < 0.89, NA, j))
        cov_mean_coef <- as.vector(direct_coef_means_nonspat[x,])
        cov_summary <- na.omit(cbind(cov_prop_retained, cov_mean_coef))
        cleaned_prop_retained <- cov_prop_retained[!is.na(cov_prop_retained)]
        cov_df <- data.frame(Proportion_retained = cleaned_prop_retained,
                             Mean_coefficient = cov_summary[, 2])
        cov_df <- cov_df[order(-abs(cov_df[, 2])),]
      })
      names(key_covariates) <- colnames(data)[1:n_nodes]
      
      #Calculate mean indirect interaction estimates from bootstrapped models
      indirect_coef_list <- booted_mrfs %>%
        purrr::map('indirect_coefs')
      
      rm(booted_mrfs)
      
      indirect_coef_means <-lapply(seq_along(indirect_coef_list[[1]]), function(x){
        Reduce(`+`, sapply(indirect_coef_list, "[[", x)) / length(indirect_coef_list)
      })
      names(indirect_coef_means) <- names(indirect_coef_list[[1]])
      
      #Clear un-needed objects and return parameters of interest
      gc()
      list(key_covariates = key_covariates,
           raw_coefs = direct_coef_list,
           indirect_coefs = indirect_coef_means)
    }, cl = cl)
    stopCluster(cl)
    
  } else {
    
    #### If parallel loading fails, or if n_cores = 1, use lapply instead ####
    cat('Fitting bootstrap_MRF models in sequence using 1 core... \n')
    
    lambda_results <- pbapply::pblapply(lambda1_seq, function(l) {
      booted_mrfs <- lapply(seq_len(n_bootstraps[[l]]), function(x) {
        sample_data <- sample(seq_len(100), 1)
        
        if(spatial){
          ID <- as.data.frame(id_data$Animal.ID)
          colnames(ID) <- 'Animal.ID'
          booted_data <- shuffle_rows(cbind(data, ID))
          sample_data <- booted_data[, 1:(ncol(booted_data) - ncol(ID))]
          
          if(nas_present){
            sample_data <- impute_nas(sample_data)
          }
          
          sample_data <- prep_MRF_covariates(sample_data, n_nodes)
          sample_id <- booted_data[, ((ncol(booted_data) + 1) -
                                            ncol(ID)):ncol(booted_data)]
          sample_id <- as.data.frame(sample_id)
          colnames(sample_id) <- 'Animal.ID'
          sample_id$Animal.ID <- as.factor(sample_id$Animal.ID)
          
          # Use invisible to prevent printing of timing messages in each iteration
          invisible(utils::capture.output(mod <- MRFcov_ID(data = sample_data,
                                                                symmetrise = symmetrise,
                                                                id_data = sample_id,
                                                                n_nodes = n_nodes,
                                                                n_cores = 1,
                                                                prep_covariates = FALSE,
                                                                n_covariates = n_covariates,
                                                                family = family,
                                                                bootstrap = TRUE)))
          
        } else {
          sample_data <- shuffle_rows(data)
          
          if(nas_present){
            sample_data <- impute_nas(booted_data)
          }
          
          sample_data <- prep_MRF_covariates(sample_data, n_nodes)
          invisible(utils::capture.output(mod <- MRFcov(data = sample_data,
                                                        symmetrise = symmetrise,
                                                        n_nodes = n_nodes,
                                                        n_cores = 1,
                                                        prep_covariates = FALSE,
                                                        n_covariates = n_covariates,
                                                        family = family,
                                                        bootstrap = TRUE)))
        }
        
        list(direct_coefs = mod$direct_coefs,
             indirect_coefs = mod$indirect_coefs)
      })
      
      direct_coef_list <- booted_mrfs %>%
        purrr::map('direct_coefs')
      
      direct_coef_means <- apply(array(unlist(direct_coef_list),
                                       c(nrow(booted_mrfs[[1]]$direct_coefs),
                                         ncol(booted_mrfs[[1]]$direct_coefs),
                                         length(booted_mrfs))), c(1, 2), mean)
      rownames(direct_coef_means) <- rownames(booted_mrfs[[1]]$direct_coefs)
      colnames(direct_coef_means) <- colnames(booted_mrfs[[1]]$direct_coefs)
      
      #Calculate proportion of bootstrap samps in which each direct coefficient occurs
      n_total_covariates <- ncol(booted_mrfs[[1]]$direct_coefs)
      
      prop_covs_retained <- matrix(0, n_nodes, n_total_covariates)
      for(i in seq_len(n_nodes)){
        for(j in seq_len(n_total_covariates)){
          prop_covs_retained[i,j] <- 1 - countzero(booted_mrfs, i, j)
        }
      }
      rownames(prop_covs_retained) <- rownames(booted_mrfs[[1]]$direct_coefs)
      colnames(prop_covs_retained) <- colnames(booted_mrfs[[1]]$direct_coefs)
      
      # Remove spatial splines from the key covariates list
      if(spatial){
        prop_covs_retained <- prop_covs_retained[, !grepl('Spatial', colnames(prop_covs_retained))]
        direct_coef_means_nonspat <- direct_coef_means[, !grepl('Spatial', colnames(direct_coef_means))]
      } else {
        direct_coef_means_nonspat <- direct_coef_means
      }
      
      #Create list of covariates that are not zero for each node in data
      key_covariates <- lapply(seq_len(n_nodes), function(x){
        cov_prop_retained <- apply(data.frame(prop_covs_retained)[x,], 2,
                                   function(j) ifelse(j < 0.9, NA, j))
        cov_mean_coef <- as.vector(direct_coef_means_nonspat[x,])
        cov_summary <- na.omit(cbind(cov_prop_retained, cov_mean_coef))
        cleaned_prop_retained <- cov_prop_retained[!is.na(cov_prop_retained)]
        cov_df <- data.frame(Proportion_retained = cleaned_prop_retained,
                             Mean_coefficient = cov_summary[, 2])
        cov_df <- cov_df[order(-abs(cov_df[, 2])), ]
      })
      names(key_covariates) <- colnames(data)[1:n_nodes]
      
      #Calculate mean indirect interaction coefficients from bootstrap samps
      indirect_coef_list <- booted_mrfs %>%
        purrr::map('indirect_coefs')
      rm(booted_mrfs)
      
      indirect_coef_means <-lapply(seq_along(indirect_coef_list[[1]]), function(x){
        Reduce(`+`, sapply(indirect_coef_list, "[[", x)) / length(indirect_coef_list)
      })
      names(indirect_coef_means) <- names(indirect_coef_list[[1]])
      rm(indirect_coef_list)
      
      #Clear un-needed objects and return parameters of interest
      gc()
      list(key_covariates = key_covariates,
           raw_coefs = direct_coef_list,
           indirect_coefs = indirect_coef_means)
    })
  }
  
  #### Calculate summary statistics of coefficients from bootstrapped models ####
  #Name each list element by its iteration value
  names(lambda_results) <- lambda1_seq
  
  #Calculate summary coefficient statistics across all iterations
  all_direct_coef_list <- lambda_results %>%
    purrr::map('raw_coefs') %>%
    purrr::flatten()
  
  all_direct_coef_means <- apply(array(unlist(all_direct_coef_list),
                                       c(nrow(lambda_results[[1]]$raw_coefs[[1]]),
                                         ncol(lambda_results[[1]]$raw_coefs[[1]]),
                                         total_reps)), c(1, 2), mean)
  rownames(all_direct_coef_means) <- rownames(lambda_results[[1]]$raw_coefs[[1]])
  colnames(all_direct_coef_means) <- colnames(lambda_results[[1]]$raw_coefs[[1]])
  
  # Remove spatial splines from the mean coefficients list for calculating relative importance
  if(spatial){
    all_direct_coef_means_nonspat <- all_direct_coef_means[, !grepl('Spatial', colnames(all_direct_coef_means))]
  } else {
    all_direct_coef_means_nonspat <- all_direct_coef_means
  }
  
  all_direct_coef_upper90 <- apply(array(unlist(all_direct_coef_list),
                                         c(nrow(lambda_results[[1]]$raw_coefs[[1]]),
                                           ncol(lambda_results[[1]]$raw_coefs[[1]]),
                                           total_reps)), c(1, 2),
                                   function(x){quantile(x, probs = 0.95)})
  rownames(all_direct_coef_upper90) <- rownames(lambda_results[[1]]$raw_coefs[[1]])
  colnames(all_direct_coef_upper90) <- colnames(lambda_results[[1]]$raw_coefs[[1]])
  
  all_direct_coef_lower90 <- apply(array(unlist(all_direct_coef_list),
                                         c(nrow(lambda_results[[1]]$raw_coefs[[1]]),
                                           ncol(lambda_results[[1]]$raw_coefs[[1]]),
                                           total_reps)), c(1, 2),
                                   function(x){quantile(x, probs = 0.05)})
  rownames(all_direct_coef_lower90) <- rownames(lambda_results[[1]]$raw_coefs[[1]])
  colnames(all_direct_coef_lower90) <- colnames(lambda_results[[1]]$raw_coefs[[1]])
  
  #Calculate relative importance of key covariates across the set of iterations
  
  coef_rel_importances <- t(apply(all_direct_coef_means_nonspat[, -1], 1, function(i) i^2 / sum(i^2)))
  mean_key_coefs <- lapply(seq_len(n_nodes), function(x){
    if(length(which(coef_rel_importances[x, ] > 0.01)) == 1){
      node_coefs <- data.frame(Variable = names(which((coef_rel_importances[x, ] > 0.01) == T)),
                               Rel_importance = coef_rel_importances[x, which(coef_rel_importances[x, ] > 0.01)],
                               Mean_coef = all_direct_coef_means_nonspat[x, 1 + which(coef_rel_importances[x, ] > 0.01)])
    } else {
      node_coefs <- data.frame(Variable = names(coef_rel_importances[x, which(coef_rel_importances[x, ] > 0.01)]),
                               Rel_importance = coef_rel_importances[x, which(coef_rel_importances[x, ] > 0.01)],
                               Mean_coef = all_direct_coef_means_nonspat[x, 1 + which(coef_rel_importances[x, ] > 0.01)])
    }
    
    rownames(node_coefs) <- NULL
    
    node_coefs <- node_coefs[order(-node_coefs[, 2]), ]
  })
  names(mean_key_coefs) <- rownames(all_direct_coef_means_nonspat)
  
  #Calculate indirect coefficient summary statistics
  all_indirect_coef_list <- lambda_results %>%
    purrr::map('indirect_coefs')
  
  all_indirect_coef_means <-lapply(seq_along(all_indirect_coef_list[[1]]), function(x){
    Reduce(`+`, sapply(all_indirect_coef_list, "[", x)) / length(all_indirect_coef_list)
  })
  names(all_indirect_coef_means) <- names(all_indirect_coef_list[[1]])
  
  # If Poisson, return scale factors for back-conversion of coefficients
  if(return_poisson){
    return(list(direct_coef_means = all_direct_coef_means,
                direct_coef_upper90 = all_direct_coef_upper90,
                direct_coef_lower90 = all_direct_coef_lower90,
                indirect_coef_mean = all_indirect_coef_means,
                mean_key_coefs = mean_key_coefs,
                mod_type = 'bootstrap_MRF',
                mod_family = 'poisson',
                poiss_sc_factors = poiss_sc_factors))
    
    # Else, return list without scale factors
  } else {
    return(list(direct_coef_means = all_direct_coef_means,
                direct_coef_upper90 = all_direct_coef_upper90,
                direct_coef_lower90 = all_direct_coef_lower90,
                indirect_coef_mean = all_indirect_coef_means,
                mean_key_coefs = mean_key_coefs,
                mod_type = 'bootstrap_MRF',
                mod_family = family))
  }
  
}