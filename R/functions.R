library(dplyr)
library(icenReg)

# Assume the scores are in {1, 2, ..., score_max}
simulate_data_dep <- function(n, max_score, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k <- max_score - 1  # k is the max improvement

  # Simulate covariates
  race <- sample(1:3, n, replace = TRUE, prob = rep(1/3, 3))
  X_num <- sample(1:3, n, replace = TRUE, prob = rep(1/3, 3))

  # Define possible values for Delta and T
  delta_vals <- -k:k
  T_vals     <- 1:k

  # Pre-allocate vectors for Delta and T.
  Delta <- numeric(n)
  T <- numeric(n)

  #Sample Delta and T conditionally on the covariates
  for (i in 1:n) {

    weights_delta <- rep(1/(2*k + 1), 2*k+1)
    Delta[i] <- sample(delta_vals, 1, prob = weights_delta)

    weights_T <- rep(1/k, k)
    T[i] <- sample(T_vals, 1, prob = weights_T)
  }

  #Compute the indicator variable (feeling significantly better)
  Indicator <- as.integer(Delta >= T) #indicator whether patient feel significantly improved

  #Combine into a data frame
  # 'race' remains as a numeric variable.
  df <- data.frame(
    Delta = Delta,
    T = T,
    Indicator = Indicator,
    race = race,
    X_num = X_num
  )
  return(df)
}

#Assume the scores are in {1, 2, ..., score_max}
transform_df <- function(df, score_max){

  #check if there is non-sensicle data
  if(any(df$Delta < 1 & df$Indicator == 1)){
    warning("Exists at least one patient with no score improvement but reports feeling significantly better")
  }
  if(any(df$Delta == (score_max-1) & df$Indicator == 0)){
    warning("Exists at least one patient with max score improvement but does not feel significantly better")
  }

  df$L <- ifelse(df$Indicator == 1, 0, df$Delta) # (0, score_max - 1] not (1, score_max - 1]
  df$R <- ifelse(df$Indicator == 1, df$Delta, score_max-1)

  return(df)
}


mean_from_surv <- function(model, cov = NULL) {

  model <- getSCurves(model, cov)

  # Extract Turnbull intervals (a k x 2 matrix)
  Tbull_ints <- model$Tbull_ints

  # Extract survival curve estimates
  S_curve <- model$S_curves[[1]] #rather strange, but this is how it goes
  k <- length(S_curve)

  # Compute probability masses from the survival drops:
  p <- numeric(k)
  p[1] <- 1 - S_curve[1]
  if (k > 1) {
    p[2:k] <- S_curve[1:(k - 1)] - S_curve[2:k]
  }

  points <- Tbull_ints[, 2]
  mean_event_time <- sum(p * points)

  return(mean_event_time)
}

MIC <- function(df, method, user_formula = NULL, cov_input = NULL) {

  if(method == "np"){
    model <- ic_np(cbind(df$L, df$R))
  } else if(method %in% c("po", "ph")){
    full_formula <- as.formula(
      paste("Surv(L, R, type = 'interval2') ~", user_formula)
    )
    model <- ic_sp(full_formula, data = df, model = method)
  } else {
    stop("Unsupported method: ", method)
  }

  # Common operations for all methods
  meanT <- mean_from_surv(model)
  meanD <- mean(df$Delta)

  result <- list(
    model = model,
    meanT = meanT, #meanT can be conditioning, depending on function argument
    meanD = meanD
  )

  return(result)
}
