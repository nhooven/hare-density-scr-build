# PROJECT: Building hare SCR model
# SCRIPT: 06 - Model selection
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 05 Jan 2026
# COMPLETED: 05 Jan 2026
# LAST MODIFIED: 05 Jan 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# Here we'll compare detection functions with a Bayes factor
# According to Dey et al. 2019, the harmonic mean estimator of the marginal likelihood
# is probably sufficient.

# We can extract posterior draws of the log-likelihood using this procedure:
# https://danielturek.github.io/public/sumPostLogDens/sumPostLogDens.html

# ______________________________________________________________________________
# 2. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(nimble)
library(coda)
library(tictoc)
library(MCMCvis)

# ______________________________________________________________________________
# 3. Read in data ----
# ______________________________________________________________________________

# S polygon
S.sf <- st_read(dsn = paste0(getwd(), "/Data for model/", lyr = "focal_S.shp"))

# traps points
traps.sf <- st_read(dsn = paste0(getwd(), "/Data for model/", lyr = "focal_traps.shp"))

# mark-recapture data
load(paste0(getwd(), "/Data for model/mr_data.RData"))

# trap operation matrix
load(paste0(getwd(), "/Data for model/trap_op.RData"))

# ______________________________________________________________________________
# 4. Separate list ----
# ______________________________________________________________________________

# capture histories
ch <- focal.mr.3[[1]]

# previous capture indicator
prev.c <- focal.mr.3[[2]]

# individual index
indivs <- focal.mr.3[[3]]

# sex indicator (0 is F)
sex <- ifelse(focal.mr.3[[5]] == "F",
              0,
              1)

# ______________________________________________________________________________
# 5. Define required quantities ----
# ______________________________________________________________________________

# number of captured individuals
n.ind = nrow(ch)

# number of traps
J = nrow(traps.sf)

# number of occasions
K = ncol(ch)

# S limits
S.xlim = as.numeric(st_bbox(S.sf)[c(1, 3)])
S.ylim = as.numeric(st_bbox(S.sf)[c(2, 4)])

# trap coordinates
x.j = st_coordinates(traps.sf)[ , 1]
y.j = st_coordinates(traps.sf)[ , 2]

# ______________________________________________________________________________
# 6. Data augmentation ----

# as a general rule of thumb, we'll start adding ~2x as many augmented to the captured 
# We will obviously have to modify this for the open model

# ______________________________________________________________________________

n.aug = n.ind * 2

M = n.ind + n.aug

# data augmentation variable
z <- c(rep(1, times = n.ind), rep(NA, times = n.aug))

# add capture histories
ch.1 <- rbind(ch, matrix(J + 1,    # index for "no captures"
                         nrow = M - n.ind, 
                         ncol = ncol(ch)))

# add to previous capture covariate
prev.c.1 <- rbind(prev.c, matrix(0,
                                 nrow = M - n.ind, 
                                 ncol = ncol(prev.c)))

# ______________________________________________________________________________
# 7. Build data lists ----
# ______________________________________________________________________________

# constants
constant.list <- list(
  
  M = M,
  J = J,
  K = K,
  S.area = st_area(S.sf) * 0.0001,   # area in ha,
  n.sex = 2
  
)

# data
data.list <- list(
  
  # individual data
  ch = ch.1,
  prev.c = prev.c.1,
  sex = c(sex, rep(NA, times = M - n.ind)),     # latent for augmented
  
  # trap operation
  trap.op = trap.op,
  
  # state space
  S.xlim = S.xlim,
  S.ylim = S.ylim,
  
  # traps
  x.j = x.j,
  y.j = y.j,
  
  # data augmentation
  z = z
  
)

# ______________________________________________________________________________
# 8. Modeling ----
# ______________________________________________________________________________
# 8a. EXPONENTIAL DETECTION ----
# ______________________________________________________________________________

model.1.code <- nimbleCode({
  
  # log density
  logDens ~ dnorm(0, 1)    ## this distribution does not matter
  
  # priors
  
  # sex specific parameters
  for (t in 1:n.sex) {
    
    alpha0[t] ~ dnorm(0, sd = 3)           # on the logit scale
    alpha2[t] ~ dnorm(0, sd = 2)           # coefficient
    sigma[t] ~ dunif(0, 60)                # this should be > 20
    
    alpha1[t] <- -1 / sigma[t]
    
  }
  
  # data augmentation
  psi ~ dunif(0, 1)
  psi.sex ~ dunif(0, 1)
  
  # loop through inndividuals M
  for (i in 1:M) {
    
    # data augmentation parameter
    z[i] ~ dbern(psi)
    
    # sex
    sex[i] ~ dbern(psi.sex)
    sex2[i] <- sex[i] + 1    # integer for indexing
    
    # activity centers
    s[i, 1] ~ dunif(S.xlim[1], S.xlim[2])
    s[i, 2] ~ dunif(S.ylim[1], S.ylim[2])
    
    # loop through traps J
    for (j in 1:J) {
      
      # distance from each AC s to each trap j
      d[i, j] <- pow(pow(s[i, 1] - x.j[j], 2) + pow(s[i, 2] - y.j[j], 2), 0.5)
      
    }
    
    # loop through occasions K
    for (k in 1:K) {
      
      # and loop through traps
      for (j in 1:J) {
        
        # linear predictor (with trap operation matrix)
        lp[i, k, j] <- exp(alpha0[sex2[i]] + alpha1[sex2[i]] * d[i, j] + alpha2[sex2[i]] * prev.c[i, k]) * z[i] * trap.op[j, k]
        
        # probability
        p[i, k, j] <- lp[i, k, j] / (1 + sum(lp[i, k, 1:J]))  # sum over all traps
        
      }
      
      # probability of not being captured as the complement of all trap-specific probs
      p[i, k, J + 1] <- 1 - sum(p[i, k, 1:J])
      
      # categorical likelihood
      ch[i, k] ~ dcat(p[i, k, 1:(J + 1)])
      
    }
    
  }
  
  # derived quantities
  N <- sum(z[1:M])
  
  D <- N / S.area
  
})

# ______________________________________________________________________________
# 8b. BIVARIATE GAUSSIAN DETECTION ----
# ______________________________________________________________________________

model.2.code <- nimbleCode({
  
  # log density
  logDens ~ dnorm(0, 1)    ## this distribution does not matter
  
  # priors
  
  # sex specific parameters
  for (t in 1:n.sex) {
    
    alpha0[t] ~ dnorm(0, sd = 3)           # on the logit scale
    alpha2[t] ~ dnorm(0, sd = 2)           # coefficient
    sigma[t] ~ dunif(0, 60)                # this should be > 20
    
    alpha1[t] <- -1 / (2 * sigma[t] * sigma[t])
    
  }
  
  # data augmentation
  psi ~ dunif(0, 1)
  psi.sex ~ dunif(0, 1)
  
  # loop through inndividuals M
  for (i in 1:M) {
    
    # data augmentation parameter
    z[i] ~ dbern(psi)
    
    # sex
    sex[i] ~ dbern(psi.sex)
    sex2[i] <- sex[i] + 1    # integer for indexing
    
    # activity centers
    s[i, 1] ~ dunif(S.xlim[1], S.xlim[2])
    s[i, 2] ~ dunif(S.ylim[1], S.ylim[2])
    
    # loop through traps J
    for (j in 1:J) {
      
      # distance from each AC s to each trap j
      d[i, j] <- pow(pow(s[i, 1] - x.j[j], 2) + pow(s[i, 2] - y.j[j], 2), 0.5)
      
    }
    
    # loop through occasions K
    for (k in 1:K) {
      
      # and loop through traps
      for (j in 1:J) {
        
        # linear predictor (with trap operation matrix)
        lp[i, k, j] <- exp(alpha0[sex2[i]] + alpha1[sex2[i]] * pow(d[i, j], 2) + alpha2[sex2[i]] * prev.c[i, k]) * z[i] * trap.op[j, k]
        
        # probability
        p[i, k, j] <- lp[i, k, j] / (1 + sum(lp[i, k, 1:J]))  # sum over all traps
        
      }
      
      # probability of not being captured as the complement of all trap-specific probs
      p[i, k, J + 1] <- 1 - sum(p[i, k, 1:J])
      
      # categorical likelihood
      ch[i, k] ~ dcat(p[i, k, 1:(J + 1)])
      
    }
    
  }
  
  # derived quantities
  N <- sum(z[1:M])
  
  D <- N / S.area
  
})

# ______________________________________________________________________________
# 8c. Initial values ----
# ______________________________________________________________________________

inits <- list(
  
  # initial s - this will just be uniform draws
  s = cbind(runif(M, S.xlim[1], S.xlim[2]),
            runif(M, S.ylim[1], S.ylim[2])),
  
  sigma = c(25, 25),
  
  alpha0 = c(-2, -2),
  
  alpha2 = c(0, 0),
  
  psi = 0.5,
  
  psi.sex = 0.5,
  
  z = c(rep(NA, times = n.ind), rep(0, times = M - n.ind)),
  
  sex = c(rep(NA, times = n.ind), rep(0, times = M - n.ind))
  
)

# ______________________________________________________________________________
# 8d. Parameters to monitor ----
# ______________________________________________________________________________

monitor <- c(
  
  "sigma",
  "alpha0",
  "alpha2",
  "psi",
  "psi.sex",
  "D",
  "s",
  "sex",
  "logDens"
  
)

# ______________________________________________________________________________
# 8e. Function to save the posterior log-density ----
# ______________________________________________________________________________

sumLogPostDens <- nimbleFunction(
  
  name = 'sumLogPostDens',
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    stochNodes <- setdiff(model$getNodeNames(stochOnly = TRUE), target)
    
  },
  
  run = function() {
    
    model[[target]] <<- model$getLogProb(stochNodes)
    
  },
  
  methods = list( reset = function() {} )
)

# ______________________________________________________________________________
# 8f. Set up and run models ----
# ______________________________________________________________________________

# EXPONENTIAL DETECTION
# set up model
model.1 <- nimbleModel(
  
  code = model.1.code,
  constants = constant.list,
  data = data.list,
  inits = inits
  
)

# add custom sampler
conf <- configureMCMC(model.1, print = FALSE)
conf$removeSamplers('logDens')   ## remove sampler assigned to 'logDens'
conf$addSampler(target = 'logDens', type = 'sumLogPostDens')   ## add our custom sampler
conf$printSamplers()

# build model
model.1.mcmc <- buildMCMC(
  
  conf = model.1,
  monitors = monitor
  
)

# compile model
compileNimble(model.1)
model.1.comp <- compileNimble(model.1.mcmc)

# run MCMC
tic()
model.1.run <- runMCMC(
  
  mcmc = model.1.comp,
  niter = 10000,
  nburnin = 2000,
  nchains = 3,
  samplesAsCodaMCMC = TRUE
  
)
toc()

# BIVARIATE GAUSSIAN DETECTION
# set up model
model.2 <- nimbleModel(
  
  code = model.2.code,
  constants = constant.list,
  data = data.list,
  inits = inits
  
)

# add custom sampler
conf <- configureMCMC(model.2, print = FALSE)
conf$removeSamplers('logDens')   ## remove sampler assigned to 'logDens'
conf$addSampler(target = 'logDens', type = 'sumLogPostDens')   ## add our custom sampler
conf$printSamplers()

# build model
model.2.mcmc <- buildMCMC(
  
  conf = model.2,
  monitors = monitor
  
)

# compile model
compileNimble(model.2)
model.2.comp <- compileNimble(model.2.mcmc)

# run MCMC
tic()
model.2.run <- runMCMC(
  
  mcmc = model.2.comp,
  niter = 10000,
  nburnin = 2000,
  nchains = 3,
  samplesAsCodaMCMC = TRUE
  
)
toc()

# ______________________________________________________________________________
# 9. Summarize and visualize ----
# ______________________________________________________________________________

# model 1
MCMCsummary(model.1.run,
            params = c("logDens"))

hist(as.data.frame(model.1.run[[1]])$logDens)

# model 2
MCMCsummary(model.2.run,
            params = c("logDens"))

hist(as.data.frame(model.2.run[[1]])$logDens)

# calculate harmonic means
hm_estimator <- function (samples) {
  
  # extract summed log densities, and exponentiate
  dens <- c(exp(as.data.frame(samples[[1]])$logDens),
            exp(as.data.frame(samples[[2]])$logDens),
            exp(as.data.frame(samples[[3]])$logDens))
  
  # harmonic mean
  hm <- 1 / ((1 / length(dens)) * sum(1 / dens))
  
  # return
  return(hm)
  
}

# calculate Bayes factor
hm_estimator(model.1.run) / hm_estimator(model.2.run)

# weak evidence that the BvG works better
MCMCsummary(model.1.run, params = c("D"))
MCMCsummary(model.2.run, params = c("D"))

# D estimates are about the same

# we'll revisit this with a larger model
