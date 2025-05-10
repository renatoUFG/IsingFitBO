# GNU GENERAL PUBLIC LICENSE - Version 2, June 1991
# Copyright (C) Claudia van Borkulo
# Modified by Renato Rodrigues Silva (2025) to include Bayesian Optimization.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
IsingFitBO = function(x, method="BayesOpt", family = "binomial", 
                      AND = TRUE, niter=20, plot = TRUE,
                      gamma_hyp = 0.25, ...)   
{
  t0 = Sys.time()
  xx = x
  set.seed(123)
  if (family != "binomial") {
    stop("This procedure is only supported for binary 
         (family='binomial') data")
  }
  log_uniform_grid = function(center, ngrid = 6, log_range = 3) {
    # Generates a grid of values regularly spaced 
    #on a logarithmic scale around a given center point. 
    #The center does not need to be included
    #in the grid, but the grid is symmetric around log(center).
    log_center = log(center)
    log_min = log_center - log_range / 2
    log_max = log_center + log_range / 2
    log_seq = seq(log_min, log_max, length.out = ngrid)
    return(exp(log_seq))
  }
  # The `allowedNodes` function checks whether a 
  # binary variable (node) has sufficient variability 
  # to be included in logistic regression modeling 
  # within the IsingFit procedure.
  allowedNodes <- function(nodeValues) {
    nodeValues <- as.factor(nodeValues)
    valuesFrequency <- table(nodeValues)
    minFrequency <- min(valuesFrequency)
    maxFrequency <- max(valuesFrequency)
    if (minFrequency <= 1 || maxFrequency >= length(nodeValues) - 1) {
      return(0)
    } else {
      return(1)
    }
  }
  #allthemeans = colMeans(x)
  NodesToAnalyze = apply(x, 2, allowedNodes) != 0
  names(NodesToAnalyze) = colnames(x)
  if (!any(NodesToAnalyze)) stop("No variance in dataset")
  if (any(!NodesToAnalyze)) {
    warning(paste("Nodes with too little variance (not allowed):", 
                  paste(colnames(x)[!NodesToAnalyze], 
                        collapse = ", ")))
  }
  x = as.matrix(x)
  # Check data: Any missing?
  if (any(is.na(x))){
    stop("IsingFitBO does not support missing data")
  }
  # Binary?
  if (!all( x == 1 | x== 0)) {
    stop("IsingFit only supports data that is encoded as {0,1}")
  }
  x = x[, NodesToAnalyze, drop = FALSE]
  nvar = ncol(x)
  nobs = nrow(x)
  p = nvar - 1
  EBIC_per_node = numeric(nvar)
  lambda_per_node = numeric(nvar)
  intercepts_per_node = vector("list", length=nvar)
  betas_per_node = vector("list", length=nvar)
  if(method=="BayesOpt"){
    for(i in 1:nvar){
      if(nobs >= p){
        Lambdas = c(0.00001,seq(0.01, by = 0.01, length=7),3)
      } else {
        Lambdas = c(0.001,seq(0.01, by = 0.01, length=7),3)
      }
      intercepts = betas =  lambdas = vector("list", length=nvar)
      sumloglik <- J <- EBIC <-  numeric(length(Lambdas))
      for (iter in 1:niter) {
        #Step 1 - Fit Poisson LASSO Regression 
        #for each ith - node and k-th lambda
        mod = glmnet(x[, -i], 
                     x[,i], lambda = Lambdas,
                     family = "binomial")
        if(iter > 1 ){
          intercepts[[i]] = c(intercepts[[i]],mod$a0) 
          betas[[i]] = cbind(betas[[i]],mod$beta)
          lambdas[[i]] =  c(lambdas[[i]],mod$lambda) 
          Lambdas = lambdas[[i]]
        } else {
          intercepts[[i]] = mod$a0
          betas[[i]] = mod$beta
          lambdas[[i]] =  mod$lambda
        }
        #Step 2 - Compute EBIC for each ith node 
        #and kth lambda
        if(iter == 1){
          y = x[,i]
          eta <- predict(mod, newx = x[, -i], 
                             s = lambdas[[i]], 
                             type = "link")
          mu = 1 / (1 + exp(-eta))
          sumloglik = colSums(y * log(mu) + (1 - y) * log(1 - mu))
          J = colSums(betas[[i]] != 0)
          EBIC = -2 * sumloglik + J * log(nobs) + 2 * gamma_hyp * J * log(p)
        } else {
          k = length(Lambdas)
          eta = predict(mod, 
                        newx = x[,-i], 
                        s = lambdas[[i]][k], 
                        type = "link")
          mu = 1 / (1 + exp(-eta))  
          y = x[,i]
          sumloglik[k] = sum(y * log(mu) + (1 - y) * log(1 - mu))
          beta_k = betas[[i]][,k]
          J[k] = sum(beta_k != 0)
          EBIC[k] = -2 * sumloglik[k] + J[k] * log(nobs) + 
            2 * gamma_hyp * J[k] * log(p)
        } 
        #Step 3 - Obtain the minimal value of EBIC
        y_min = min(EBIC)
        lambda_best = lambdas[[i]][which.min( EBIC)]
        #Step 4 - Generate new values of lambda
        center <- lambda_best
        Lambdas_prime = log_uniform_grid(center, 
                                         ngrid = length(lambdas[[i]]), 
                                         log_range = 2)
        #Step 5 - Fit GP for regression
        gp_model =  gpkm(lambdas[[i]], EBIC, 
                         kernel = "Gaussian",
                         nug.est = FALSE, nug=1e-06)
        #Step 6 - Estimate the Expected Improvement
        pred <- gp_model$pred(Lambdas_prime, se.fit = TRUE)
        mu <- pred$mean
        sigma <- pred$se
        Z = (y_min - mu) / sigma
        EI = (y_min - mu) * pnorm(Z) + sigma * dnorm(Z)
        lambda_EI = Lambdas_prime[which.max(EI)]
        Lambdas = lambda_EI
      } 
      EBIC_per_node[[i]] = min(EBIC)
      lambda_per_node[[i]] = lambdas[[i]][which.min(EBIC)]
      intercepts_per_node[[i]] = intercepts[[i]][which.min(EBIC)]
      betas_per_node[[i]] = betas[[i]][,which.min(EBIC)]
    }  
  } else if(method == "Grid"){
    for(i in 1:nvar){
      if(nobs >= p){
        Lambdas = c(0.00001,seq(0.01, by = 0.01, length=7),3)
      } else {
        Lambdas = c(0.001,seq(0.01, by = 0.01, length=7),3)
      }
      intercepts = betas =  lambdas = vector("list", length=nvar)
      sumloglik <- J <- EBIC <-  numeric(length(Lambdas))
      #Step 1 - Fit Poisson LASSO Regression 
      #for each ith - node and k-th lambda
      mod = glmnet(x[, -i], 
                   x[,i], lambda = Lambdas,
                   family = "binomial")
      intercepts[[i]] = mod$a0
      betas[[i]] = mod$beta
      lambdas[[i]] =  mod$lambda
      #Step 2 - Compute EBIC for each ith node 
      #and kth lambda
      y = x[,i]
      eta_mat <- predict(mod, newx = x[, -i], 
                         s = lambdas[[i]], 
                         type = "link")
      mu_mat <- 1 / (1 + exp(-eta_mat))
      logliks <- colSums(y * log(mu_mat) + (1 - y) * log(1 - mu_mat))
      J <- colSums(mod$beta != 0)
      EBIC <- -2 * logliks + J * log(nobs) + 2 * gamma_hyp * J * log(p)
      #Step 3 - Choose the value of lambda that minimizes
      #EBIC
      EBIC_per_node[[i]] = min(EBIC)
      lambda_per_node[[i]] = lambdas[[i]][which.min(EBIC)]
      intercepts_per_node[[i]] = intercepts[[i]][which.min(EBIC)]
      betas_per_node[[i]] = betas[[i]][,which.min(EBIC)]
    }    
  } else {
    stop("The function only accepts the strings 'BayesOpt' 
         or 'Grid' in the method argument.")
  } 
  EBIC = EBIC_per_node
  lambdas = lambda_per_node
  betas = betas_per_node
  thresholds = intercepts_per_node
  weights.opt <- matrix(, nvar, nvar)
  for (i in 1:nvar) {
    weights.opt[i, -i] = betas[[i]]
  }
  asymm.weights = weights.opt
  diag(asymm.weights) = 0
  if (AND == TRUE) {
    adj = weights.opt
    adj = (adj != 0) * 1
    EN.weights = adj * t(adj)
    EN.weights = EN.weights * weights.opt
    meanweights.opt = (EN.weights + t(EN.weights)) / 2
    diag(meanweights.opt) = 0
  } else {
    meanweights.opt = (weights.opt + t(weights.opt)) / 2
    diag(meanweights.opt) = 0
  }
  graphNew <- matrix(0, length(NodesToAnalyze), length(NodesToAnalyze))
  graphNew[NodesToAnalyze, NodesToAnalyze] <- meanweights.opt
  colnames(graphNew) <- rownames(graphNew) <- colnames(xx)
  threshNew = rep(NA, p)
  threshNew[NodesToAnalyze] = thresholds
  if (plot == TRUE) notplot <- FALSE else notplot <- TRUE
  q = qgraph(graphNew, layout = "spring", 
              labels = names(NodesToAnalyze), 
              DoNotPlot = notplot)
  Res <- list(
    weiadj = graphNew, 
    thresholds = threshNew, 
    q = q, gamma = gamma_hyp,
    AND = AND, 
    time = Sys.time() - t0, 
    asymm.weights = asymm.weights,
    lambda.values = lambdas,
    EBIC = EBIC
  )
  class(Res) <- "IsingFitBO"
  return(Res)
}   
     





