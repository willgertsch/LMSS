# generate log-likelihood function for use with optimizer
#
ll_factory = function(type, Y, X, Z) {

  dimZ = ncol(Z)
  dimX = ncol(X)
  if (type == 'continuous') {
    ll_fun = function(var) {
      # beta0 = var[1:2]
      # beta1 = var[3:4]
      # gamma = var[5:(4+ncol(X))]
      # sigma = var[(5+ncol(X)):length(var)]

      beta0 = var[1:dimZ]
      beta1 = var[(dimZ + 1):(2*dimZ)]
      gamma = var[(2*dimZ + 1):(2*dimZ + dimX)]
      sigma = var[(2*dimZ + dimX + 1):(2*dimZ + dimX + 2)]

      # check constraint
      # looking for treatment effect difference
      if (beta1[2] <= 0)
        return(-Inf)

      eta = X %*% gamma
      p = ilogit(eta)
      mu0 = Z %*% beta0
      mu1 = Z %*% (beta0 + beta1)
      f0 = dnorm(Y, mu0, sigma[1])
      f1 = dnorm(Y, mu1, sigma[2])
      ll = sum(log(p * f1 + (1 - p) * f0))
      if (is.na(ll))
        return(-Inf)
      else
        return(ll)
    }
  }
  else if (type == 'binary') {

    ll_fun = function(var) {

      # beta0 = var[1:2]
      # beta1 = var[3:4]
      # gamma = var[5:(4+ncol(X))]
      beta0 = var[1:dimZ]
      beta1 = var[(dimZ + 1):(2*dimZ)]
      gamma = var[(2*dimZ + 1):(2*dimZ + dimX)]

      # check constraint
      if (beta1[2] <= 0)
        return(-Inf)

      eta = X %*% gamma
      p = ilogit(eta)
      mu0 = Z %*% beta0
      mu1 = Z %*% (beta0 + beta1)
      f0 = dbinom(Y, 1, ilogit(mu0))
      f1 = dbinom(Y, 1, ilogit(mu1))
      ll = sum(log(p * f1 + (1 - p) * f0))
      if (is.na(ll))
        return(-Inf)
      else
        return(ll)
    }
  }
  else if (type == 'count') {
    ll_fun = function(var) {

      # beta0 = var[1:2]
      # beta1 = var[3:4]
      # gamma = var[5:(4+ncol(X))]
      beta0 = var[1:dimZ]
      beta1 = var[(dimZ + 1):(2*dimZ)]
      gamma = var[(2*dimZ + 1):(2*dimZ + dimX)]

      # check constraint
      if (beta1[2] <= 0)
        return(-Inf)

      eta = X %*% gamma
      p = ilogit(eta)
      mu0 = Z %*% beta0
      mu1 = Z %*% (beta0 + beta1)
      f0 = dpois(Y, exp(mu0))
      f1 = dpois(Y, exp(mu1))
      ll = sum(log(p * f1 + (1 - p) * f0))
      if (is.na(ll))
        return(-Inf)
      else
        return(ll)
    }
  }
  else
    stop('ll_factory: type not supported')

  return(ll_fun)

}
