# logistic normal mixture model

ilogit = function(eta) {
  1/(1 + exp(-eta))
}


# log-likelihood
lnm_logl = function(beta1, beta2, sigma1, sigma2, gamma, Y, Z, X) {


  eta = X %*% gamma
  p = ilogit(eta)

  mu1 = Z %*% (beta1 + beta2)
  mu2 = Z %*% beta1

  f1 = dnorm(Y, mu1, sigma1)
  f2 = dnorm(Y, mu2, sigma2)

  ll = sum(log(p * f1 + (1 - p) * f2))

  return(ll)
}

# generate a sample from the logistic normal mixture model given mean parameters
# n: sample size
# mu1: Z * (beta1 + beta2)
# mu2: Z * beta1
# eta: X * gamma
# Return: outcome vector Y
rlnm = function(n, mu1, mu2, eta, sigma1, sigma2) {

  # randomly generate binary class based on logistic regression
  class = rbinom(n, 1, ilogit(eta))

  # sample from mixture model with known class
  Y = numeric(n)
  for (i in 1:n) {
    if (class[i] == 1) {
      Y[i] = rnorm(1, mu1[i], sigma1)
    }
    else if (class[i] == 0) {
      Y[i] = rnorm(1, mu2[i], sigma2)
    }
  }

  # return outcome
  return(Y)

}

# update the value of gamma during M step using NR
# X: matrix of covariates
# W: vector of posterior probabilities of class membership
# gamma_init: initial gamma values
# iter: number of iterations for NR
# Returns: updated gamma
max_Q1 = function(X, W, gamma_init, iter, silent) {

  gamma = gamma_init # initial value
  for (i in 1:iter) {
    P = ilogit(X %*% gamma)
    grad = t(X) %*% (W - P)
    Omega = diag(c(P * (1 - P)))
    IM = t(X) %*% Omega %*% X
    gamma = gamma + solve(IM) %*% grad
    if (!silent)
      cat("Iteration: ", i,  ",  Gamma=", gamma, "\n")
  }
  return(gamma)
}

# update beta parameters and sigma using WLS steps
# Y: outcome vector
# Z: matrix of covariates
# W: Vector of posterior probabilities of class membership
# *_init: initial values for parameters
# Return: updated beta1, beta2, sigma
max_Q2 = function(Y, Z, W, beta1_init, beta2_init, sigma1_init, sigma2_init) {

  tbeta1_init = beta1_init + beta2_init # reparameterize
  tbeta2_init = beta1_init

  # beta estimation using WLS-like step
  # update tbeta1 and then tbeta2
  f1 = exp(-1/(2*sigma1_init^2) * (Y - Z %*% tbeta1_init)^2)
  f2 = exp(-1/(2*sigma2_init^2) * (Y - Z %*% tbeta2_init)^2)
  omega1 = 1/(W * f1 + (1 - W) * f2) * f1
  Omega1 = diag(c(omega1))
  hat_tbeta1 = solve(t(Z) %*% Omega1 %*% Z) %*% t(Z) %*% Omega1 %*% Y
  sigma2_1 = sum(omega1*(Y - Z%*%hat_tbeta1)^2)/sum(omega1)
  hat_sigma1 = sqrt(sigma2_1)

  #f1 = exp(-1/(2*sigma1_init^2) * (Y - Z %*% hat_tbeta1)^2)
  f1 = exp(-1/(2*hat_sigma1^2) * (Y - Z %*% hat_tbeta1)^2)

  omega2 = 1/(W * f1 + (1 - W) * f2) * f2
  Omega2 = diag(c(omega2))
  hat_tbeta2 = solve(t(Z) %*% Omega2 %*% Z) %*% t(Z) %*% Omega2 %*% Y

  # sigma estimation using WLS-like step
  # shortcut: averaging individual variance estimates
  #sigma2_1 = sum(omega1*(Y - Z%*%hat_tbeta1)^2)/sum(omega1)
  sigma2_2 = sum(omega2*(Y - Z%*%hat_tbeta2)^2)/sum(omega2)
  #hat_sigma1 = sqrt(sigma2_1)
  hat_sigma2 = sqrt(sigma2_2)

  return(list(beta1 = hat_tbeta2, beta2 = hat_tbeta1 - hat_tbeta2, sigma1 = hat_sigma1, sigma2 = hat_sigma2))

}

# perform one step of the EM update
# returns list of updated parameter values
update_EM = function(Y, Z, X, beta1, beta2, sigma1, sigma2, gamma, silent) {

  eta = X %*% gamma
  p = ilogit(eta)

  # E step
  mu1 = Z %*% (beta1 + beta2)
  mu2 = Z %*% beta1
  f1 = dnorm(Y, mu1, sigma1)
  f2 = dnorm(Y, mu2, sigma2)
  denoms = p*f1 + (1-p)*f2
  W = p * f1 / denoms

  # M step ####
  # gamma update
  # how many iterations of IRLS is good?
  gamma_t = max_Q1(X, W, gamma, 10, silent)
  # update beta and sigma
  result = max_Q2(Y, Z, W, beta1, beta2, sigma1, sigma2)

  # return update values
  return(list(
    beta1 = result$beta1,
    beta2 = result$beta2,
    gamma = gamma_t,
    sigma1 = result$sigma1,
    sigma2 = result$sigma2
  ))

}

# fit the logistic normal model using metaheuristics
# Y: outcome vector
# Z: design matrix for mixture regression
# X: design matrix for logistic regression
# iter: maximum number of iterations for the algorithm
# swarm: population size for algorithm
# algorithm: one of the options from metaheuristicOpt
# Return: a named list of parameters
# This function assumes that the second column of Z contains the treatment variable
# therefore, we make the assumption that beta2[2] > 0 for identifibility
# also assume that there intercepts in both regression equations
fit_lnm = function(Y, Z, X, iter, swarm, algorithm) {

  # get dimensions
  n = length(Y)
  q1 = ncol(Z)
  q2 = ncol(X)

  # set up objective function
  obj_fun = function(param) {

    beta1 = param[1:q1]
    beta2 = param[(q1 + 1):(2*q1)]
    gamma = param[(2*q1+1):(2*q1 + q2)]
    sigma1 = param[2*q1 + q2 + 1]
    sigma2 = param[2*q1 + q2 + 2]
    LL = lnm_logl(beta1, beta2, sigma1, sigma2, gamma, Y, Z, X)

    # deal with missing
    if (is.na(LL))
      return(-Inf)
    else
      return(LL)
  }

  # set bounds for problem
  minY = min(Y)
  maxY = max(Y)
  rangeY = maxY - minY
  beta1_bounds = matrix(c(
    minY, maxY, rep(c(-rangeY, rangeY), q1 - 1)
  ), nrow = 2)
  beta2_bounds = matrix(c(
    -rangeY, rangeY, 0, rangeY, rep(c(-rangeY, rangeY), q1 - 2)
  ), nrow = 2)
  gamma_bounds = matrix(c(
    -10, 10, rep(c(-5, 5), q2 - 1)
  ), nrow = 2)
  sigma_bounds = matrix(c(
    0.1, sd(Y), 0.1, sd(Y)
  ), nrow = 2)

  bounds = cbind(beta1_bounds, beta2_bounds, gamma_bounds, sigma_bounds)

  # call metaheuristics library to maximize likelihood
  control = list()
  control$maxIter = iter
  control$numPopulation = swarm
  out = metaheuristicOpt::metaOpt(
    FUN = obj_fun,
    optimType = "MAX",
    algorithm = algorithm,
    numVar = 2*q1 + q2 + 2,
    rangeVar = bounds,
    control = control
  )


  # extract and return parameters
  beta1 = out$result[1:q1]
  beta2 = out$result[(q1 + 1):(2*q1)]
  gamma = out$result[(2*q1+1):(2*q1 + q2)]
  sigma = out$result[(2*q1 + q2 + 1):(2*q1 + q2 + 2)]
  return(list(
    ll = out$optimumValue[1],
    beta1 = beta1,
    beta2 = beta2,
    gamma = gamma,
    sigma1 = sigma[1],
    sigma2 = sigma[2]
  ))

}

# EM algorithm implementation
lnm_EM = function(Y, Z, X, maxiter, silent = F) {

  # get dimensions
  n = length(Y)
  q1 = ncol(Z)
  q2 = ncol(X)

  # generate initial parameter values
  # some random
  beta1 = solve(t(Z) %*% Z) %*% t(Z) %*% Y # LS estimate
  beta2 = rnorm(q1, 0, sd(Y))
  beta2[2] = abs(beta2[2]) # beta22 > 0
  gamma = rnorm(q2, 0, 2) # values from roughly +-5
  sigma1 = sd(Y)
  sigma2 = sd(Y)

  # initial update
  result = update_EM(Y, Z, X, beta1, beta2, sigma1, sigma2, gamma, silent)
  obj = lnm_logl(result$beta1, result$beta2, result$sigma1, result$sigma2, result$gamma, Y, Z, X)

  obj_old = obj
  beta1 = result$beta1
  beta2 = result$beta2
  gamma = result$gamma
  sigma1 = result$sigma1
  sigma2 = result$sigma2


  # main loop
  for (iter in 1:maxiter) {



    # update
    result = update_EM(Y, Z, X, beta1, beta2, sigma1, sigma2, gamma, silent)
    obj = lnm_logl(result$beta1, result$beta2, result$sigma1, result$sigma2, result$gamma, Y, Z, X)

    # print iteration number of objective value
    if (!silent)
      cat("Iter: ", iter, " obj: ", obj, "\n")

    # check convergence criterion
    # will also stop if negative increase in ll
    ftolrel = 1e-12
    if (obj - obj_old < ftolrel*(abs(obj_old)+1)) {
      # case when EM step decreased ll
      if (obj - obj_old < 0)
        break
      # otherwise accept results
      obj_old = obj
      beta1 = result$beta1
      beta2 = result$beta2
      gamma = result$gamma
      sigma1 = result$sigma1
      sigma2 = result$sigma2
      break
    }
    obj_old = obj
    beta1 = result$beta1
    beta2 = result$beta2
    gamma = result$gamma
    sigma1 = result$sigma1
    sigma2 = result$sigma2

  }

  # return
  return(list(
    ll = obj_old,
    beta1 = beta1,
    beta2 = beta2,
    gamma = gamma,
    sigma1 = sigma1,
    sigma2 = sigma2
  ))
}

# predict which latent class an individual is assigned to
# X: design matrix of individual characteristics
# gamma: gamma estimate from lnm model
# Return: n x 2 matrix of class predictions and probabilities
predict_class = function(X, gamma) {

  eta = X %*% gamma
  p = ilogit(eta)
  class = ifelse(p > 0.5, 1, 0)
  return(cbind(class, p))

}
