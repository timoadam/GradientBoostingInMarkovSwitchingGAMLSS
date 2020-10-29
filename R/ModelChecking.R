## Compute the log-forward probabilities
lforward = function(mod) {
  mus = cbind(as.vector(fitted(mod$mod[[1]])$mu), as.vector(fitted(mod$mod[[2]])$mu))
  sigmas = cbind(as.vector(fitted(mod$mod[[1]])$sigma), as.vector(fitted(mod$mod[[2]])$sigma))
  n = length(mod$y)
  lalpha = matrix(NA, 2, n)
  P = dnorm(mod$y[1], mean = mus[1,], sd = exp(sigmas[1,]))
  foo = mod$delta * P
  sumfoo = sum(foo)
  lscale = log(sumfoo)
  foo = foo / sumfoo
  lalpha[, 1] = lscale + log(foo)
  for(i in 2:n) {
    P = dnorm(mod$y[i], mean=mus[i,], sd = exp(sigmas[i,]))
    foo = foo %*% mod$gamma * P
    sumfoo = sum(foo)
    lscale = lscale + log(sumfoo)
    foo = foo / sumfoo
    lalpha[, i] = log(foo) + lscale
  }
  return(lalpha)
}

## Compute the one-step-ahead forecast pseudo-residuals
PseudoResiduals <- function(mod) {
  mus = cbind(as.vector(fitted(mod$mod[[1]])$mu), as.vector(fitted(mod$mod[[2]])$mu))
  sigmas = cbind(as.vector(fitted(mod$mod[[1]])$sigma), as.vector(fitted(mod$mod[[2]])$sigma))
  la = t(lforward(mod = mod))
  n = length(mod$y)
  Res = rep(NA, n)
  pMat = matrix(NA, nrow = n, ncol = 2)
  pMat[1,] = pnorm(mod$y[1], mean=mus[1,], sd = exp(sigmas[1,]))
  Res[1] = qnorm(t(mod$delta) %*% pMat[1,])
  for(i in 2:n) {
    pMat[i,] = pnorm(mod$y[i], mean=mus[i,], sd=exp(sigmas[i,]))
    c = max(la[i - 1,])
    a = exp(la[i - 1,] - c)
    Res[i] = qnorm(t(a) %*% (mod$gamma / sum(a)) %*% pMat[i,])
  }
  return(list(Res = Res))
}
