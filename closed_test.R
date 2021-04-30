#######################################
# closed testing procedure to select validation model
# source: Vergouwe Y, Nieboer D, Oostenbrink R, et al. A closed testing procedure to select an appropriate method for updating prediction models. Stat Med 2017; 36: 4529–39.
ClosedTest <- function(coefs, X, y){
  # Implement closed testing procedure (Version: 11-01-2013)
  # Arguments:
  # coefs: Vector containing the regression coefï¬cients of the model that
  # is updated.
  # X: predictor matrix
  # y: outcome vector
  # Results:
  # coef_new: regression coefï¬cients of chosen model
  require(rms)
  if(class(X)=="data.frame"){
    X <- data.matrix(X)
  }
  if(ncol(X)!=(length(coefs)-1)){
    stop("Number of predictors not equal to the number of coefï¬cients")
  }
  n_coefs <- length(coefs)
  lp_old <- X %*% as.matrix(coefs[2:n_coefs])
  # Calculate updated model intercept
  intercept <- lrm.fit(y = y, offset = lp_old)$coefficients
  coefs_int <- c(intercept , coefs[2:n_coefs])
  # Calculate coefï¬cients after recalibration
  recal <- lrm.fit(x = lp_old, y = y)$coefficients
  coefs_recal <- c(recal[1], recal[2] * coefs[2:n_coefs])
  # Calculate coefï¬cients after model revision
  coefs_refit <- lrm.fit(x = X, y = y)$coefficients
  # Calculate the log-likelihood of the different models
  lp <- cbind(1, X) %*% coefs
  ll_original <- sum(y * lp - log(1 + exp(lp)))
  lp <- cbind(1, X) %*% coefs_int
  ll_intercept <- sum(y * lp - log(1 + exp(lp)))
  lp <- cbind(1, X) %*% coefs_recal
  ll_recalibration <- sum(y * lp - log(1 + exp(lp)))
  lp <- cbind(1, X) %*% coefs_refit
  ll_revision <- sum(y * lp - log(1 + exp(lp)))
  # Calculate difference in log-likelihood for testing of the models
  dev_original <- -2 * ll_original + 2 * ll_revision
  dev_intercept <- -2 * ll_intercept + 2 * ll_revision
  dev_recalibration <- -2 * ll_recalibration + 2 * ll_revision
  # See if difference in model ï¬t was signiï¬cant
  test1 <- (1 - pchisq(dev_original, ncol(X) + 1)) < 0.05
  test2 <- (1 - pchisq(dev_intercept, ncol(X))) < 0.05
  test3 <- (1 - pchisq(dev_recalibration, ncol(X) - 1)) < 0.05
  # See which model is chosen, index_test indicates the chosen model
  # 1. Original model
  # 2. Model with updated intercept
  # 3. Recalibrated model
  # 4. Revised model
  test_original <- 1 * (!test1)
  test_intercept <- 2 * ((test1)&(!test2))
  test_recalibration <- 3 * ((test1)&(test2)&(!test3))
  test_revision <- 4 * ((test1)&(test2)&(test3))
  index_test <- (test_original + test_intercept + test_recalibration +test_revision)
  coefs_result <- rbind(coefs, coefs_int, coefs_recal, coefs_refit)
  # Output of the function
  new_coefs <- coefs_result[index_test, ]
  model <- c("Original Model", "Model with updated intercept",
             "Recalibrated model", "Model Revision")[index_test]
  cat("Method chosen by closed test procedure:\n", model, "\n",
      "Resulting coefï¬cients:\n", new_coefs, "\n")
  res <- list(model = model, coefs = coefs_result)
  return(res)
}