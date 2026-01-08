# PROJECT: Building hare SCR model
# SCRIPT: 08 - 2022f model
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Jan 2026
# COMPLETED: 
# LAST MODIFIED: 08 Jan 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# we'll use the Richard Chandler trick here

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

# S polygons
S.list <- list(st_read(dsn = paste0(getwd(), "/Data for model/S/", lyr = "S_2A.shp")),
               st_read(dsn = paste0(getwd(), "/Data for model/S/", lyr = "S_2B.shp")),
               st_read(dsn = paste0(getwd(), "/Data for model/S/", lyr = "S_2C.shp")),
               st_read(dsn = paste0(getwd(), "/Data for model/S/", lyr = "S_3A.shp")),
               st_read(dsn = paste0(getwd(), "/Data for model/S/", lyr = "S_3C.shp")))

# traps points
traps.list <- list(st_read(dsn = paste0(getwd(), "/Data for model/traps/", lyr = "traps_2A.shp")),
                   st_read(dsn = paste0(getwd(), "/Data for model/traps/", lyr = "traps_2B.shp")),
                   st_read(dsn = paste0(getwd(), "/Data for model/traps/", lyr = "traps_2C.shp")),
                   st_read(dsn = paste0(getwd(), "/Data for model/traps/", lyr = "traps_3A.shp")),
                   st_read(dsn = paste0(getwd(), "/Data for model/traps/", lyr = "traps_3C.shp")))

# mark-recapture data
load(paste0(getwd(), "/Data for model/MR/mr_data_2A.RData"))
load(paste0(getwd(), "/Data for model/MR/mr_data_2B.RData"))
load(paste0(getwd(), "/Data for model/MR/mr_data_2C.RData"))
load(paste0(getwd(), "/Data for model/MR/mr_data_3A.RData"))
load(paste0(getwd(), "/Data for model/MR/mr_data_3C.RData"))

# trap operation matrix
load(paste0(getwd(), "/Data for model/Trap op/trap_op_2A.RData"))
load(paste0(getwd(), "/Data for model/Trap op/trap_op_2B.RData"))
load(paste0(getwd(), "/Data for model/Trap op/trap_op_2C.RData"))
load(paste0(getwd(), "/Data for model/Trap op/trap_op_3A.RData"))
load(paste0(getwd(), "/Data for model/Trap op/trap_op_3C.RData"))

# ______________________________________________________________________________
# 4. Separate lists ----
# ______________________________________________________________________________
# 4a. Bind in capture histories ----

# capture histories n x K
# we need to add NA columns for additional occasions (max = 10)

# function
# pass in a list of mr data
bind_ch <- function (mr.list) {
  
  max.K = max(unlist(lapply(mr.list, ncol)))
  
  ch <- list()
  
  for (i in 1:length(mr.list)) {
    
    ch[[i]] <- cbind(mr.list[[i]],
                      matrix(data = NA,
                             nrow = nrow(mr.list[[i]]),
                             ncol = max.K - ncol(mr.list[[i]])))
    
  }
  
  ch.mat <- do.call(rbind, ch)
  
  return(ch.mat)
  
}

# ______________________________________________________________________________

ch <- bind_ch(mr.list = list(focal.mr.2A.2[[1]],
                             focal.mr.2B.2[[1]],
                             focal.mr.2C.2[[1]],
                             focal.mr.3A.2[[1]],
                             focal.mr.3C.2[[1]]))

# ______________________________________________________________________________
# 4b. Previous capture indicator ----

# we'll use the same function and different inputs

# ______________________________________________________________________________

prev.c <- bind_ch(mr.list = list(focal.mr.2A.2[[2]],
                                 focal.mr.2B.2[[2]],
                                 focal.mr.2C.2[[2]],
                                 focal.mr.3A.2[[2]],
                                 focal.mr.3C.2[[2]]))

# ______________________________________________________________________________
# 4c. Individual indices ----
# ______________________________________________________________________________

# we;ll have to cross-reference this with the open CHs later
indivs <- 1:nrow(ch)

# group index 
group <- c(rep(1, times = length(focal.mr.2A.2[[3]])),
           rep(2, times = length(focal.mr.2B.2[[3]])),
           rep(3, times = length(focal.mr.2C.2[[3]])),
           rep(4, times = length(focal.mr.3A.2[[3]])),
           rep(5, times = length(focal.mr.3C.2[[3]])))


# sex indicator (0 is F)
sex <- ifelse(c(focal.mr.2A.2[[5]],
                focal.mr.2B.2[[5]],
                focal.mr.2C.2[[5]],
                focal.mr.3A.2[[5]],
                focal.mr.3C.2[[5]]) == "F",
              0,
              1)

# ______________________________________________________________________________
# 5. Define required quantities ----

# quite a bit of this can be packed into a separate script later
# the S bounds, trap locations, etc. can all be pre-attributed
# especially since they only vary by g and not g x t

# ______________________________________________________________________________
# 5a. Constants ----
# ______________________________________________________________________________

# number of sessions
G = 5

# number of captured individuals (by session)
n.ind = c(nrow(focal.mr.2A.2[[1]]),
          nrow(focal.mr.2B.2[[1]]),
          nrow(focal.mr.2C.2[[1]]),
          nrow(focal.mr.3A.2[[1]]),
          nrow(focal.mr.3C.2[[1]]))

# number of traps (always 36)
J = 36

# number of occasions (by session!)
K = c(ncol(focal.mr.2A.2[[1]]),
      ncol(focal.mr.2B.2[[1]]),
      ncol(focal.mr.2C.2[[1]]),
      ncol(focal.mr.3A.2[[1]]),
      ncol(focal.mr.3C.2[[1]]))

# ______________________________________________________________________________
# 5b. S limits by site ----

# function
S_limits <- function(S.list) {
  
  # define blank matrices - each [G x 2]
  S.xlim <- matrix(data = NA, nrow = G, ncol = 2)
  S.ylim <- matrix(data = NA, nrow = G, ncol = 2)
  
  # loop through
  for (i in 1:length(S.list)) {
    
    S.xlim[i, ] = as.numeric(st_bbox(S.list[[i]])[c(1, 3)])
    S.ylim[i, ] = as.numeric(st_bbox(S.list[[i]])[c(2, 4)])
    
  }
  
  return(list(S.xlim, S.ylim))
  
}

# ______________________________________________________________________________

S.limits <- S_limits(S.list)

# ______________________________________________________________________________
# 5c. Trap coordinates by site ----

# function
j_coords <- function(traps.list) {
  
  # trap coordinates - [J x 2 x G]
  j.coords <- array(data = NA, dim = c(J, 2, G))
  
  # loop through
  for (i in 1:G) {
    
    j.coords[ , 1, i] = st_coordinates(traps.list[[i]])[ , 1]
    j.coords[ , 2, i] = st_coordinates(traps.list[[i]])[ , 2]
    
  }
  
  return(j.coords)
  
}

# ______________________________________________________________________________

j.coords <- j_coords(traps.list)

# ______________________________________________________________________________
# 6. Data augmentation ----
# ______________________________________________________________________________

n.aug <- vector(length = G) 

n.aug[1] = n.ind[1] * 4
n.aug[2] = n.ind[2] * 4
n.aug[3] = n.ind[3] * 4
n.aug[4] = n.ind[4] * 4
n.aug[5] = n.ind[5] * 4

M = sum(n.ind) + sum(n.aug)

# data augmentation variable
z <- c(rep(1, times = sum(n.ind)), rep(NA, times = sum(n.aug)))

# zeroes variable
zeroes <- c(rep(NA, times = sum(n.ind)), rep(0, times = sum(n.aug)))

# add to previous capture covariate
prev.c.1 <- rbind(prev.c, matrix(0,
                                 nrow = M - sum(n.ind), 
                                 ncol = ncol(prev.c)))

# the zeroes in the 9 and 10 columns are probably okay here since we're only 
# looping through K for each session (MIGHT HAVE TO CHANGE)

# group assignment
group.1 <- c(group,
             rep(1, times = n.aug[1]),
             rep(2, times = n.aug[2]),
             rep(3, times = n.aug[3]),
             rep(4, times = n.aug[4]),
             rep(5, times = n.aug[5]))

# sex
sex.1 <- c(sex, rep(NA, times = M - sum(n.ind)))

# ______________________________________________________________________________
# 7. Build data lists ----
# ______________________________________________________________________________

# constants
constant.list <- list(
  
  # scalar constants
  n0 = sum(n.ind),
  M = M,
  J = J,
  G = G,
  n.sex = 2,
  
  # group-specific constants
  S.area = c(st_area(S.list[[1]]) * 0.0001,   # area in ha,
             st_area(S.list[[2]]) * 0.0001,
             st_area(S.list[[3]]) * 0.0001,
             st_area(S.list[[4]]) * 0.0001,
             st_area(S.list[[5]]) * 0.0001),
  
  # known indices
  K = K,
  group = group.1
  
)

# data
data.list <- list(
  
  # individual data [n0 x max(K)]
  ch = ch,
  
  # individual data [M x max(K)]
  prev.c = prev.c.1,
  sex = sex.1,
  
  # trap operation (WILL NEED TO MODIFY TO MAKE MATRICES CONFORM TO MAX(K))
  # okay to hard code for now
  trap.op = array(c(trap.op.2A, 
                    trap.op.2B,
                    trap.op.2C,
                    trap.op.3A,
                    trap.op.3C), dim = c(J, 8, G)),
  
  # state space [G x 2]
  S.xlim = S.limits[[1]],
  S.ylim = S.limits[[2]],
  
  # traps
  j.coords = j.coords,
  
  # data augmentation [M]
  z = z,
  zeroes = zeroes,
  
  # which groups to turn z on and off for N sum
  which.groups.1 = ifelse(group.1 == 1, 1, 0),
  which.groups.2 = ifelse(group.1 == 2, 1, 0),
  which.groups.3 = ifelse(group.1 == 3, 1, 0),
  which.groups.4 = ifelse(group.1 == 4, 1, 0),
  which.groups.5 = ifelse(group.1 == 5, 1, 0)
  
)

# ______________________________________________________________________________
# 8. Modeling ----
# ______________________________________________________________________________
# 8a. Code ----
# ______________________________________________________________________________

model.1.code <- nimbleCode({
  
  # priors
  # sex specific parameters [n.sex]
  for (t in 1:n.sex) {
    
    alpha0[t] ~ dnorm(0, sd = 2)           # on the logit scale
    alpha2[t] ~ dnorm(0, sd = 2)           # coefficient
    sigma[t] ~ dunif(10, 60)               # this should be > 20
    
    alpha1[t] <- -1 / sigma[t]
    
  }
  
  # data augmentation [G]
  for (g in 1:G) {
    
    psi[g] ~ dunif(0, 1)
    psi.sex[g] ~ dunif(0, 1)
    
  }
  
  # loop through all individuals [M]
  for (i in 1:M) {
    
    # data augmentation parameter (by group)
    z[i] ~ dbern(psi[group[i]])
    
    # sex
    sex[i] ~ dbern(psi.sex[group[i]]) # (by group)
    sex2[i] <- sex[i] + 1    # integer for indexing
    
    # activity centers [M, 2] (index by group)
    s[i, 1] ~ dunif(S.xlim[group[i] , 1], S.xlim[group[i] , 2])
    s[i, 2] ~ dunif(S.ylim[group[i] , 1], S.ylim[group[i] , 2])
    
    # loop through traps J
    for (j in 1:J) {
      
      # distance from each AC s to each trap j [M, J] (index by group, dim 3)
      d[i, j] <- sqrt(pow(s[i, 1] - j.coords[j, 1, group[i]], 2) + pow(s[i, 2] - j.coords[j, 2, group[i]], 2))
      
    }
    
    # loop through occasions K (index by group)
    for (k in 1:(K[group[i]])) {
      
      # and loop through traps J
      for (j in 1:J) {
        
        # linear predictor 
        lp[i, k, j] <- exp(
          
          # baseline detection by sex
          alpha0[sex2[i]] + 
            
            # distance term
            alpha1[sex2[i]] * d[i, j] + 
            
            # previous capture covariate
            alpha2[sex2[i]] * prev.c[i, k]) * 
          
          # data augmentation
          z[i] *
          
          # trap operation [J, K[g], g]
          trap.op[j, k, group[i]]
        
        # probability
        p[i, k, j] <- lp[i, k, j] / (1 + sum(lp[i, k, 1:J]))  # sum over all traps
        
      }
      
      # probability of not being captured as the complement of all trap-specific probs
      p[i, k, J + 1] <- 1 - sum(p[i, k, 1:J])
      
    }
    
  }
  
  # likelihood for the observed individuals
  for (i in 1:n0) {
    
    # loop through occasions K (index by group)
    for (k in 1:(K[group[i]])) {
      
      ch[i, k] ~ dcat(p[i, k, 1:(J + 1)])
      
    }
    
  }
  
  # likelihood for the non-detected individuals
  # inside, the probability of at AT LEAST ONE detection
  for (i in (n0 + 1):M) {
    
    zeroes[i] ~ dbern(1 - prod(1 - p[i, 1:(K[group[i]]), 1:J]))
    
  }
  
  # derived quantities [G]
  # N individuals
  N[1] <- sum(z[1:M] * which.groups.1[1:M])
  N[2] <- sum(z[1:M] * which.groups.2[1:M])
  N[3] <- sum(z[1:M] * which.groups.3[1:M])
  N[4] <- sum(z[1:M] * which.groups.4[1:M])
  N[5] <- sum(z[1:M] * which.groups.5[1:M])
  
  # density D
  D[1:G] <- N[1:G] / S.area[1:G]
  
})

# ______________________________________________________________________________
# 8b. Initial values ----
# ______________________________________________________________________________

inits <- list(
  
  # initial s - start at 0 ,0
  s = cbind(runif(constant.list$M, -100, 100),
            runif(constant.list$M, -100, 100)),
  
  alpha0 = runif(2, -5, 5),
  
  alpha2 = rnorm(2, 0, 2),
  
  sigma = runif(2, 10, 60),
  
  psi = runif(constant.list$G, 0, 1),
  
  psi.sex = runif(constant.list$G, 0, 1),
  
  z = c(rep(NA, times = sum(n.ind)), rep(0, times = M - sum(n.ind))),
  
  sex = c(rep(NA, times = sum(n.ind)), rep(0, times = M - sum(n.ind)))
  
)

# ______________________________________________________________________________
# 8c. Parameters to monitor ----
# ______________________________________________________________________________

monitor <- c(
  
  "sigma",
  "alpha0",
  "alpha2",
  "psi",
  "D"
  
)

# ______________________________________________________________________________
# 8d. Set up and run model - exponential detection ----
# ______________________________________________________________________________

# set up model
model.1 <- nimbleModel(
  
  code = model.1.code,
  constants = constant.list,
  data = data.list,
  inits = inits
  
)

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
  nburnin = 5000,
  nchains = 1,
  thin = 5,
  samplesAsCodaMCMC = TRUE
  
)
toc()

# ______________________________________________________________________________
# 9. Summarize and visualize ----
# ______________________________________________________________________________

# summary of main parameters
MCMCsummary(model.1.run,
            params = c("psi"))

# traceplot
MCMCtrace(model.1.run,
          params = "D",
          pdf = F)

# bad chain behavior, but the model is doing what I want it to do

# again, the state space is probably too small