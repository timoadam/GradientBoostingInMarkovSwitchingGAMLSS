## Loading required packages
# install.packages("fda")
# install.packages("boot")
# install.packages("numDeriv")
# install.packages("gamlss")
# install.packages("gamboostLSS")
library(fda)
library(boot)
library(numDeriv)
library(gamlss)
library(gamboostLSS)

## Function that fits Markov-switching GAMLSS using the msgamboostLSS algorithm
FitMarkovSwitchingGAMLSS <- function(x, y, N = 2, type = "MSGAMLSS", stat = FALSE, m.stop = rep(200, 2), max.iter = 50, conv.tol = 1e-03, conv.print = TRUE) { 
  x1 = x[, 1]
  delta = NULL
  gamma = matrix(c(0.95, 0.05, 0.05, 0.95), ncol = 2)
  mod = fv.mu = fv.sigma = list()
  term = FALSE
  old = 0
  circ.crit = 0
  for(j in 1:N) {
    fv.mu[[j]] = rep(c(3, 7)[j], length(y))
    fv.sigma[[j]] = rep(3, length(y))
  }
  while(term == FALSE) {
    for(i in 1:max.iter) {
      delta.next = delta
      gamma.next = gamma
      fv.mu.next = fv.mu
      fv.sigma.next = fv.sigma
      if(is.null(delta)) {
        delta = solve(t(diag(N) - gamma + 1), rep(1, N))
      }
      allprobs = matrix(NA, length(y), N)
      lalpha = lbeta = matrix(NA, N, length(y))
      for(j in 1:N) {
        allprobs[,j] = apply(as.matrix(y), 2, dNO, mu = fv.mu.next[[j]], sigma=exp(fv.sigma.next[[j]]))
      }
      allprobs = ifelse(!is.na(allprobs), allprobs, 1)
      foo = delta * allprobs[1,]
      sumfoo = sum(foo)
      lscale = log(sumfoo)
      foo = foo / sumfoo
      lalpha[, 1] = log(foo) + lscale
      for(j in 2:length(y)) {
        foo = foo %*% gamma * allprobs[j,]
        sumfoo = sum(foo)
        lscale = lscale + log(sumfoo)
        foo = foo / sumfoo
        lalpha[, j] = log(foo) + lscale
      }
      foo = rep(1 / N, N)
      lbeta[, length(y)] = rep(0, N)
      lscale = log(N)
      foo = foo / sum(foo)
      for(j in (length(y) - 1):1) {
        foo = gamma %*% (allprobs[j + 1,] * foo)
        lbeta[, j] = log(foo) + lscale
        sumfoo = sum(foo)
        foo = foo / sumfoo
        lscale = lscale + log(sumfoo)
      }                   
      lallprobs = log(allprobs)
      llh = max(lalpha[, length(y)]) + log(sum(exp(lalpha[, length(y)] - max(lalpha[, length(y)]))))
      weights = matrix(NA, N, length(y))
      for(j in 1:length(y)) {
        weights[,j] = exp(lalpha[,j] + lbeta[, j] - llh)
      }
      for(j in 1:N) {
        for(k in 1:N) {
          gamma.next[j, k] = gamma[j, k] * sum(exp(lalpha[j, 1:(length(y) - 1)] + lallprobs[2:length(y), k] + lbeta[k, 2:length(y)] - llh))
        }
      }
      gamma.next = gamma.next / apply(gamma.next, 1, sum)
      if(stat == TRUE) {
        delta.next = solve(t(diag(N) - gamma.next + 1), rep(1, N))
      }else{
        delta.next = exp(lalpha[, 1] + lbeta[, 1] - llh)
        delta.next = delta.next / sum(delta.next)
      }
      conv.crit = sum(abs(gamma - gamma.next)) + sum(abs(delta - delta.next))
      ind = weights
      for(j in 1:N){
        if(type == "MSGLMLSS") {
          mod[[j]] = glmboostLSS(y ~ x1, weights = ind[j,], data = data.frame(x1 = x1, y = y), families = as.families("NO"), method = "noncyclic", control = boost_control(mstop = m.stop[j], nu = 0.1))
        }
        if(type == "MSGAMLSS") {
          mod[[j]] = gamboostLSS(y ~ x1, weights = ind[j,], data = data.frame(x1 = x1, y = y), families = as.families("NO"), method = "noncyclic", control = boost_control(mstop = m.stop[j], nu = 0.1))
        }
        fv.mu.next[[j]] = as.vector(fitted(mod[[j]]$mu))
        fv.sigma.next[[j]] = as.vector(fitted(mod[[j]]$sigma))
      }
      conv.crit = abs(llh - old)
      if(conv.print == TRUE) {
        cat("Iteration = ", i, ", criterion = ", round(conv.crit, 3), "\r")
      }
      if(conv.crit < conv.tol | (conv.crit < 1 & abs(conv.crit - circ.crit) < 1e-09) | i == max.iter) {
        if(i == max.iter) {
          print(paste("No convergence within", max.iter, "iterations"))
        }else{
          print(paste("Convergence after", i, "iterations"))
        }
        term = TRUE
        break
      }  
      delta = delta.next
      gamma = gamma.next
      fv.mu = fv.mu.next
      fv.sigma = fv.sigma.next
      old = llh
    }
  }
  return(list(x1 = x1, y = y, delta = delta.next, gamma = gamma.next, mod = mod, m.stop = m.stop, llh = llh, state.probs = weights))
}
