# Require packages
require(metafor)
require(userfriendlyscience)

# Enter data manually
dat <- data.frame(
  yi = c(0.25, 0.61, 0.29, 0.69, 0.61, -0.31, 0.44),
  n1 = c(155, 44, 305, 59, 401, 27, 318),
  n2 = c(168, 52, 292, 62, 398, 28, 61),
  ci_lower = c(0.03, 0.2, 0.13, 0.32, 0.47, -0.84, 0.16)
)

dat$ci.halfwidth = dat$yi - dat$ci_lower
dat$sei = dat$ci.halfwidth/1.96

dat2 <- escalc(yi = yi, sei = sei, data = dat, measure = "SMD")

### fit fixed-effects model
res <- rma(yi, vi, data=dat2, measure="SMD")
forest(res)

# Model bias
# The following functions are from https://github.com/nicebread/meta-showdown

estimate.onestep.selection.heterogeneous <- function(z.obs,n1,n2,alpha,theta.init){
  p.value <- 1 - pnorm(z.obs)
  sel <- z.obs > 0 & p.value < alpha
  
  if(sum(sel)==0 | sum(sel)==length(sel)){
    if(sum(sel)==length(sel)){ w <- 1/(length(sel)+2) }		# Not identified; use Wilson-like
    if(sum(sel)==0){ w <- 1 - 1/(length(sel)+2) }			# estimator.
    e <- estimate.known.onestep.selection.heterogeneous(z.obs,n1,n2,w,alpha, theta.init[1:2])
    tmp <- e[[2]]
    tmp <- cbind(tmp, c(NA,NA))
    tmp <- rbind(tmp, c(NA,NA,NA))
    return(list(est=c(e[[1]],w), est.var=tmp, ll=e[[3]]))
  }
  
  if(sum(sel) > 0 & sum(sel)<length(sel)){
    theta.init <- c(theta.init[1], log(theta.init[2:3]))
    tmpf <- function(theta,z.obs,n1,n2,alpha){ onestep.heterogeneous.nll(c(theta[1],exp(theta[2:3])),z.obs,n1,n2,alpha) }
    tmpg <- function(theta,z.obs,n1,n2,alpha){ onestep.heterogeneous.nll(c(theta[1],theta[2:3]),z.obs,n1,n2,alpha) }
    tmpo <- optim(theta.init, tmpf, z.obs=z.obs, n1=n1, n2=n2, alpha=alpha)
    theta.hat <- c(tmpo$par[1], exp(tmpo$par[2:3]))
    tmpv <- matrix(NA, 3, 3)
    suppressWarnings(try( tmpv <- solve(optimHess(theta.hat, tmpg, z.obs=z.obs, n1=n1, n2=n2, alpha=alpha)), silent=TRUE))
    return(list(est=theta.hat, est.var=tmpv, ll=-tmpo$value))
  }	
}

# ---------------------------------------------------------------------
#  Estimate the three-parameter selection model (3PSM) implemented in McShane et al. 2016
# For estimation functions, see 7b-selection.meta.functions.R

TPSM.est <- function(t, n1, n2, long=TRUE) {	
  
  # Three-parameter selection model
  # McShane et al. implmentation: init.value gives an initial guess for the effect size, heterogeneity, and relative
  # likelihood of reporting a study that is not statistically significant and directionally consistent.
  
  # use very general starting values for the parameters (not tuned to our specific simulations)
  theta.init <- c(expected.d = 0.3, max.tau= 0.5, p.report = 0.99)
  alpha <- 0.05
  
  mm <- NULL
  try(mm <- estimate.onestep.selection.heterogeneous(t, n1, n2, alpha/2, theta.init), silent=FALSE)
  
  # initialize as empty df
  res.wide <- data.frame(
    method = "3PSM",
    term = "b0",
    estimate = NA,
    std.error = NA,
    statistic = NA,
    p.value = NA,
    conf.low = NA,
    conf.high = NA
  )
  if(!is.null(mm)) {
    if (all(complete.cases(mm[[2]]))) {
      SE <- sqrt(diag(mm[[2]]))
      res.wide <- data.frame(
        method = "3PSM",
        term = c("b0", "max.tau", "p.report"),
        estimate = mm[[1]],
        std.error = SE,
        statistic = NA,
        p.value = pnorm(abs(mm[[1]]) / SE, lower.tail=FALSE)*2,
        conf.low = mm[[1]] + qnorm(alpha/2)*SE,
        conf.high = mm[[1]] + qnorm(1-alpha/2)*SE
      )
    }
  }
  
  return(res.wide)
}

# set.seed(5)
# dat <- dataMA(k = 12, delta = 0.41, tau = 0.01, empN = TRUE, maxN=500, minN=0, meanN=0, selProp = 0, qrpEnv = "none")
# dat <- data.frame(dat)

onestep.heterogeneous.nll <- function(theta, z.obs, n1, n2, alpha=0.025){
  delta <- theta[1]						# True population average effect size
  tau <- theta[2]							# Heterogeneity
  w <- theta[3]							# Relative likelihood non-stat sig study is reported
  z.cut <- qnorm(1-alpha)
  d.obs <- z.obs / sqrt((n1*n2)/(n1+n2))
  d.cut <- z.cut / sqrt((n1*n2)/(n1+n2))
  k <- sum(z.obs<z.cut)
  
  s <- sqrt(tau^2 + 1/n1 + 1/n2)
  ll <- ifelse(k==0, 0, k*log(w))			# Equivalent to defining 0*log(0)=0 as is common
  ll <- ll + sum(dnorm(d.obs, delta, s, log=TRUE))
  ll <- ll - sum(log( w*pnorm(d.cut,delta,s) + (1 - pnorm(d.cut,delta,s)) ))
  -ll
}

### End of code from Schönbrodt


# Convert effect sizes to t before analysis
dat$t <- convert.d.to.t(d = dat$yi, n1 = dat$n1, n2 = dat$n2)

# Initialise values for analysis
theta.init <- c(expected.d = 0.41, max.tau= 0.5, p.report = 0.90)
alpha <- 0.05

# Perform analysis
est <- TPSM.est(t=dat$t, n1=dat$n1, n2=dat$n2)
est

# draw funnel plots
funnel(res, main = "Reported")
funnel(res, level=95, refline=0, main = "What I suspect")
funnel(res, level=95, refline=est$estimate[1], main = "3-process bias model") # Estimate from below

