# model fitting function
# type: continuous, binary, count
# Y: outcome
# X: subgroup design matrix => treatment indicator should always be second
# Z: outcome model design matrix
# swarm: population size for algorithm
# maxIter: maximum number of iterations
# algorithm: algorithm to use with metaOpt
# rangeVar: parameter bounds, can manually set otherwise default bounds will be set
fit_LMM = function(type, Y, X, Z, swarm = 40, maxIter = 500, algorithm, rangeVar=NULL) {

  # get likelihood function
  ll_fun = ll_factory(type, Y, X, Z)

  dimZ = ncol(Z)
  dimX = ncol(X)

  if (dimZ > 3 | dimZ < 2)
    stop("Z model is not supported.")

  # set variable bounds
  # setting reasonable bounds for logistic variables.
  # basing this on Bayesian priors for logistic regression
  # http://www.stat.columbia.edu/~gelman/research/published/priors11.pdf
  logOR_lb = -10
  logOR_ub = 5 # more than a 0.01 to 0.5 jump, which is not going to happen for this data
  if (is.null(rangeVar)) {
    if (type == 'continuous') {
      minY = min(Y)
      maxY = max(Y)
      rangeY = maxY - minY
      sdY = sd(Y)

      # switch depending on ANOVA or ANCOVA
      if (dimZ == 2) {
        rangeVar = matrix(c(
          minY, -rangeY, -rangeY, 0, rep(logOR_lb, ncol(X)), 0.1, 0.1,
          maxY, rangeY, rangeY, rangeY, rep(logOR_ub, ncol(X)), sdY, sdY
        ),
        nrow = 2, byrow = T)
      }
      else if (dimZ == 3) {
        rangeVar = matrix(c(
          minY, -rangeY, -rangeY, minY, 0, -rangeY, rep(logOR_lb, ncol(X)), 0.1, 0.1,
          maxY, rangeY, rangeY, minY, rangeY, rangeY, rep(logOR_ub, ncol(X)), sdY, sdY
        ),
        nrow = 2, byrow = T)
      }

    }
    else if (type == 'binary') {

      # switch depending on ANOVA or ANCOVA
      if (dimZ == 2) {
        rangeVar = matrix(c(
          logOR_lb, logOR_lb, logOR_lb, 0, rep(logOR_lb, ncol(X)),
          logOR_ub, logOR_ub, logOR_ub, logOR_ub, rep(logOR_ub, ncol(X))
        ),
        nrow = 2, byrow = T)
      }
      else if (dimZ == 3) {
        rangeVar = matrix(c(
          logOR_lb, logOR_lb, logOR_lb, logOR_lb, 0, logOR_lb, rep(logOR_lb, ncol(X)),
          logOR_ub, logOR_ub, logOR_ub, logOR_ub, logOR_ub, logOR_ub, rep(logOR_ub, ncol(X))
        ),
        nrow = 2, byrow = T)
      }

    }
    else if (type == 'count') {

      # think the OR bounds will also be find for Poisson regression
      # switch depending on ANOVA or ANCOVA
      if (dimZ == 2) {
        rangeVar = matrix(c(
          logOR_lb, logOR_lb, logOR_lb, 0, rep(logOR_lb, ncol(X)),
          logOR_ub, logOR_ub, logOR_ub, logOR_ub, rep(logOR_ub, ncol(X))
        ),
        nrow = 2, byrow = T)
      }
      else if (dimZ == 3) {
        rangeVar = matrix(c(
          logOR_lb, logOR_lb, logOR_lb, logOR_lb, 0, logOR_lb, rep(logOR_lb, ncol(X)),
          logOR_ub, logOR_ub, logOR_ub, logOR_ub, logOR_ub, logOR_ub, rep(logOR_ub, ncol(X))
        ),
        nrow = 2, byrow = T)
      }
    }
  }


  result = metaheuristicOpt::metaOpt(
    ll_fun,
    optimType = 'MAX',
    algorithm = algorithm,
    numVar = ncol(rangeVar),
    rangeVar = rangeVar,
    control = list(
      numPopulation = swarm,
      maxIter = maxIter
    )
  )

  # extract results
  param = as.numeric(result$result)
  # beta1 = param[1:2]
  # beta2 = param[3:4]
  # gamma = param[5:(5 + ncol(X) - 1)]
  beta0 = param[1:dimZ]
  beta1 = param[(dimZ + 1):(2*dimZ)]
  gamma = param[(2*dimZ + 1):(2*dimZ + dimX)]
  sigma = param[(2*dimZ + dimX + 1):(2*dimZ + dimX + 2)]

  if (type == 'continuous')
    sigma = tail(param, 2)
  else
    sigma = NULL


  return(list(
    beta0 = beta0,
    beta1 = beta1,
    gamma = gamma,
    sigma = sigma,
    algorithm = algorithm,
    swarm = swarm,
    maxIter = maxIter,
    ll = as.numeric(result$optimumValue),
    timeElapsed = result$timeElapsed,
    Y = Y,
    X = X,
    Z = Z,
    type = type
  ))

}
