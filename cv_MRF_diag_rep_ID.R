cv_MRF_diag_rep_ID = function(data, symmetrise, n_nodes, n_cores, id_data,
                                   sample_seed, n_folds, n_fold_runs, n_covariates,
                                   compare_null, family, plot = TRUE){
  
  #### Specify default parameter values and initiate warnings ####
  if(!(family %in% c('gaussian', 'poisson', 'binomial')))
    stop('Please select one of the three family options:
         "gaussian", "poisson", "binomial"')
  
  if(missing(symmetrise)){
    symmetrise <- 'mean'
  }
  
  if(any(!is.finite(as.matrix(data)))){
    stop('No infinite values permitted in data', call. = FALSE)
  }
  
  if(missing(compare_null)) {
    compare_null <- FALSE
  }
  
  if(missing(n_folds)) {
    if(nrow(data) < 50){
      n_folds <- 2
      warning('nrow(data) is less than 50, using 2-fold validation by default')
    } else {
      if(nrow(data) < 100){
        n_folds <- 5
        warning('nrow(data) is less than 100, using 5-fold validation by default')
      } else{
        n_folds <- 10
        warning('n_folds missing, using 10-fold validation by default')
      }
    }
  } else {
    if(sign(n_folds) == 1){
      #Make sure n_folds is a positive integer
      n_folds <- ceiling(n_folds)
    } else {
      stop('Please provide a positive integer for n_folds')
    }
  }
  
  if(missing(n_fold_runs)) {
    n_fold_runs <- n_folds
  } else {
    if(sign(n_fold_runs) == 1){
      #Make sure n_fold_runs is a positive integer
      n_fold_runs <- ceiling(n_fold_runs)
    } else {
      stop('Please provide a positive integer for n_fold_runs')
    }
  }
  
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
  
  if(missing(sample_seed)) {
    sample_seed <- ceiling(runif(1, 0, 100000))
  }
  
  #### Generate cached model(s) to avoid unneccessary refit in each run of n_fold_runs ####
  cat("Generating node-optimised spatial Conditional Random Fields model", "\n", sep = "")
  if(family == 'binomial'){
    mrf <- suppressWarnings(MRFcov_ID(data = data, 
                                      symmetrise = symmetrise, 
                                      n_nodes = n_nodes,
                                      n_cores = n_cores,
                                      family = 'binomial', 
                                      id_data = id_data))
    
    if(compare_null){
      cat("\nGenerating non-spatial model", "\n", sep = "")
      mrf_null <- suppressWarnings(MRFcov(data = data,
                                          symmetrise =  symmetrise,
                                          n_nodes = n_nodes,
                                          n_cores = n_cores,
                                          family = 'binomial'))
    }
  }
  
  if(family == 'poisson'){
    mrf <- suppressWarnings(MRFcov_ID(data = data, 
                                      symmetrise = symmetrise, 
                                      n_nodes = n_nodes,
                                      n_cores = n_cores,
                                      family = 'poisson', 
                                      id_data = id_data))
    
    if(compare_null){
      cat("\nGenerating non-spatial model", "\n", sep = "")
      mrf_null <- suppressWarnings(MRFcov(data = data,
                                          symmetrise =  symmetrise,
                                          n_nodes = n_nodes,
                                          n_cores = n_cores,
                                          family = 'poisson'))
    }
  }
  
  if(family == 'gaussian'){
    mrf <- suppressWarnings(MRFcov_spatial(data = data, 
                                           symmetrise = symmetrise, 
                                           n_nodes = n_nodes,
                                           n_cores = n_cores,
                                           family = 'gaussian', 
                                           id_data = id_data))
    
    if(compare_null){
      cat("\nGenerating non-spatial model", "\n", sep = "")
      mrf_null <- suppressWarnings(MRFcov(data = data,
                                          symmetrise =  symmetrise,
                                          n_nodes = n_nodes,
                                          n_cores = n_cores,
                                          family = 'gaussian'))
    }
  }
  
  # Store cached model(s) in a list
  cached_model <- list(mrf = mrf)
  if(compare_null){
    cached_model$mrf_null <- mrf_null
  }
  
  # Store cached predictions in a list
  cat("\nCalculating model predictions of the supplied data", "\n", sep = "")
  cat("Generating spatial MRF predictions ...", "\n", sep = "")
  cached_predictions <- list(predictions = predict_MRF(data = mrf$mrf_data,
                                                       prep_covariates = FALSE,
                                                       cached_model$mrf))
  if(compare_null){
    cat("Generating null MRF predictions ...", "\n", sep = "")
    cached_predictions$null_predictions <- predict_MRF(data, cached_model$mrf_null)
  }
  
  #### Replicate cv_MRF_diag n_fold_runs times, using the cached models in each run ####
  cat("\nCalculating predictive performance across test folds", "\n", sep = "")
  repped_cvs <- lapply(seq_len(n_fold_runs), function(x){
    cat("Processing cross-validation run ", x, " of ", n_fold_runs, " ...\n", sep = "")
    cv_MRF_diag(data = data, n_nodes = n_nodes,
                n_folds = n_folds,
                n_cores = n_cores, family = family,
                compare_null = compare_null, plot = FALSE,
                cached_model = cached_model,
                cached_predictions = cached_predictions,
                sample_seed = sample_seed,
                mod_labels = c('Spatial MRF', 'Non-spatial MRF'))
  })
  
  plot_dat <- do.call(rbind, repped_cvs)
  
  #### Return either a plot or a dataframe of predictive metrics ####
  if(plot){
    if(family == 'gaussian'){
      plot_gauss_cv_diag_optim(plot_dat, compare_null = compare_null)
    }
    
    if(family == 'poisson'){
      plot_poiss_cv_diag_optim(plot_dat, compare_null = compare_null)
    }
    
    if(family == 'binomial'){
      plot_binom_cv_diag_optim(plot_dat, compare_null = compare_null)
    }
  } else {
    return(plot_dat)
  }
}