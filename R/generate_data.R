# generate data
# N: sample size
# beta0: main effects
# beta1: subgroup interaction effects
# beta can have 2 or 3 parameters
# gamma: subgroup membership covariates
# type: outcome type, "continuous", "binary", 'count"
# sigma: for continuous outcome, set the two variance parameters
# trt_prob: treatment assignment probability
# X: subgroup covariate data.frame turned into matrix by model.matrix
generate_data = function(N, beta0, beta1, gamma, type, sigma = NULL, trt_prob, X) {


  # sanitity checks
  if (length(beta0) != length(beta1))
    stop('beta0, beta1 dimension mismatch')


  baseline_var = rnorm(N)
  # design matrix
  # match Z to dimension of beta
  if (length(beta0) == 2) {
    Z = cbind(rep(1, N), rbinom(N, 1, trt_prob))
  }
  else if (length(beta0) == 3) {
    Z = cbind(rep(1, N), rbinom(N, 1, trt_prob), baseline_var)
    X = cbind(X, baseline_var)
  }
  else
    stop("generate_data: length(beta) > 3 adding more covariates is better in real life, but for simplicity we don't.")


  # compute means
  mu0 = Z %*% beta0
  mu1 = Z %*% (beta0 + beta1)


  eta = X %*% gamma
  class = rbinom(N, 1, ilogit(eta))

  Y = numeric(N)
  if (type == 'continuous') {
    for (i in 1:N) {
      if (class[i] == 0) {
        Y[i] = rnorm(1, mu0[i], sigma[1])
      }
      else if (class[i] == 1) {
        Y[i] = rnorm(1, mu1[i], sigma[2])
      }
    }
  }
  else if (type == 'binary') {
    for (i in 1:N) {
      if (class[i] == 0) {
        Y[i] = rbinom(1, 1, ilogit(mu0[i]))
      }
      else if (class[i] == 1) {
        Y[i] = rbinom(1, 1, ilogit(mu1[i]))
      }
    }
  }
  else if (type == 'count') {

    for (i in 1:N) {
      if (class[i] == 0) {
        Y[i] = rpois(1, exp(mu0[i]))
      }
      else if (class[i] == 1) {
        Y[i] = rpois(1, exp(mu1[i]))
      }
    }
  }

  return(
    list(
      Y=Y, Z=Z, X=X, class=class
    )
  )

}

# generate subgroup covariate data frame
# converts automatically into matrix form
# N: sample size
generate_X = function(N) {
  X1 = rbinom(N, 1, 0.5)
  X2 = rnorm(N)

  X = data.frame(X1, X2)

  temp = rep(1, nrow(X))
  X = model.matrix(temp ~., data = cbind(temp, X))

  return(X)
}
