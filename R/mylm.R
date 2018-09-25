

if(FALSE){ #Part3c
  model_age <- mylm(wages ~ age, data = SLID)
  model_edu <- mylm(wages ~ education, data = SLID)
  model_both <- mylm(wages ~ education + age, data = SLID)
  model_age$coefficients
  model_edu$coefficients
  model_both$coefficients
}


if(FALSE) { #Part4
  model4a <- mylm(formula = wages ~ sex + age + language + I(education^2), data=SLID)
  summary(model4a)
  #negative intercept
  #aggregate response in education
  plot(model4a)
  model4b <- mylm(formula = wages ~ language*age, data=SLID)
  summary(model4b)
  plot(model4b)
  model4c <- mylm(formula = wages ~ -1 + education, data=SLID)
  summary(model4c)
  plot(model4c)
  #altough some of these plots do not look horrible, they all have a clear trend.
  #They should be completeley random, which they are not.
  #The intercept are clearly different from 0, indicating it should be included.
}

# Select Build, Build and reload to build and lode into the R-session.

mylm <- function(formula, data = list(), contrasts = NULL, ...){
  # Extract model matrix & responses
  mf <- model.frame(formula = formula, data = data)
  X  <- model.matrix(attr(mf, "terms"), data = mf, contrasts.arg = contrasts)
  y  <- model.response(mf)
  terms <- attr(mf, "terms")


  # Add code here to calculate coefficients, residuals, fitted values, etc...
  # and store the results in the list est
  inv_XtX <- solve(t(X) %*% X)

  # This is our estimation for \beta
  beta <- inv_XtX %*% t(X) %*% y

  # The fitted values are easily calculated
  fitted_values <- X %*% beta

  # These are the predicted values for the simplest possible model,
  # i.e. no covariates at all, only an intercept.
  Y_hat_H0 <- rep(1,length(fitted_values))*mean(y)
  residuals <- y - fitted_values
  SSE <- sum( (y - fitted_values)^2 )
  SST <- sum( (y - Y_hat_H0)^2 )
  sigma_tilde <- SSE/length(y)
  covariance_matrix <- sigma_tilde * solve(t(X) %*% X)
  #z-values and p-values
  t_value <- vector("numeric", length(beta))
  p_value <- vector("numeric", length(beta))
  for(i in 1:length(beta)){
    t_value[i] <- (beta[i] - 0)/sqrt(diag(covariance_matrix)[i])
    p_value[i] <- 2*pnorm(abs(t_value[i]), lower.tail=FALSE)

  }
  F_obs = ((length(fitted_values) - length(beta))*(SST - SSE)) / (SSE*(length(beta)-1))
  Chi <- (length(beta)-1) * F_obs
  p_chi <- pchisq(Chi, df=1, lower.tail=FALSE)
  Rsq <- 1 - SSE/SST

  est <- list(coefficients = beta, pvalues = p_value,
              zvalues = t_value, covmatrix =covariance_matrix,
              fitted_values = fitted_values, residuals = residuals,
              fstatistic = F_obs, SSE = SSE, SST = SST,
              chistatistic = Chi, p_chi = p_chi, Rsq = Rsq, model = mf)


  # Store call and formula used
  est$call <- match.call()
  est$formula <- formula

  # Set class name. This is very important!
  class(est) <- 'mylm'

  # Return the object with all results
  return(est)
}

print.mylm <- function(object, ...){
  # Code here is used when print(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Call: \n')
  print(object$call)
  cat('\nCoefficients: \n')
  temp <- matrix(c(object$coefficients[1], object$coefficients[2]), nrow=1, ncol=2)
  print_m <- data.frame(temp)
  colnames(print_m) <- c('(Intercept)', 'education')
  print_m
}

summary.mylm <- function(object, ...){
  # Code here is used when summary(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Call: \n')
  print(object$call)
  cat('\nCoefficients: \n')
  temp_matrix <- matrix( c(object$coefficients,
                           sqrt(diag(object$covmatrix)),object$zvalues,
                                object$pvalues), nrow=length(object$coefficients))
  coeff_matrix <- data.frame(temp_matrix)
  rownames(coeff_matrix) <- rownames(object$coefficients)
  colnames(coeff_matrix) <- c('Estimate', 'Std. Error', 'z-value', 'p-value')
  print(coeff_matrix)
  cat('F-statistic: ')
  print(object$fstatistic)
  cat('Chi-statistic: ')
  print(object$chistatistic)
  cat('p-value chi-statistic: ')
  print(object$p_chi)
  cat('R-squared:')
  print(object$Rsq)
}

plot.mylm <- function(object, ...){
  # Code here is used when plot(object) is used on objects of class "mylm"
  plot(object$fitted_values, object$residuals, xlab="Fitted values", ylab="Residuals")
  abline(a=0, b=0, col="gray", lty=3)
}



# This part is optional! You do not have to implement anova
anova.mylm <- function(object, ...){
  # Code here is used when anova(object) is used on objects of class "mylm"

  # Components to test
  comp <- attr(object$terms, "term.labels")

  # Name of response
  response <- deparse(object$terms[[2]])

  # Fit the sequence of models
  txtFormula <- paste(response, "~", sep = "")
  model <- list()
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula <- paste(txtFormula, comp[numComp])
    }
    else{
      txtFormula <- paste(txtFormula, comp[numComp], sep = "+")
    }
    formula <- formula(txtFormula)
    model[[numComp]] <- lm(formula = formula, data = object$model)
  }

  # Print Analysis of Variance Table
  cat('Analysis of Variance Table\n')
  cat(c('Response: ', response, '\n'), sep = '')
  cat('          Df  Sum sq X2 value Pr(>X2)\n')
  for(numComp in 1:length(comp)){
    # Add code to print the line for each model tested
  }

  return(model)

}


