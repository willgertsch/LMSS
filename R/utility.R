ilogit = function(eta) {
  1/(1 + exp(-eta))
}

# predict subgroup class based on model fit
predict_class = function(mod) {

  eta = mod$X %*% mod$gamma
  p = ilogit(eta)
  class = ifelse(p > 0.5, 1, 0)
  return(data.frame(class, p))

}

# nicely summarizes result from model
# mod: list output from fit_model
summary.subgroup = function(mod) {

  cat(mod$type, "logistic-mixture model fit using", mod$algorithm, '\n')
  cat("Log-likelihood:", mod$ll, '\n')
  cat("Estimated treatment effect for s=0:", mod$beta0[2], '\n')
  cat("Estimated treatment effect for s=1:", mod$beta0[2] + mod$beta1[2], '\n')
  cat('Treatment effect increase:', mod$beta1[2], '\n')

  class_dat = predict_class(mod)
  cat('Subgroup:', table(class_dat$class), '\n')

}
