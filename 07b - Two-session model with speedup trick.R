# PROJECT: Building hare SCR model
# SCRIPT: 07b - Two-session model with speedup trick
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Jan 2026
# COMPLETED: 06 Jan 2026
# LAST MODIFIED: 07 Jan 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# we'll use the Richard Chandler trick:
# https://groups.google.com/g/hmecology/c/o6cWDqHHgOE

# in which we add a "zeroes" variable for all augmented individuals
# and a everDetected variable for the conditional statement

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
S.2A <- st_read(dsn = paste0(getwd(), "/Data for model/S/", lyr = "S_2A.shp"))
S.2B <- st_read(dsn = paste0(getwd(), "/Data for model/S/", lyr = "S_2B.shp"))

# traps points
traps.2A <- st_read(dsn = paste0(getwd(), "/Data for model/traps/", lyr = "traps_2A.shp"))
traps.2B <- st_read(dsn = paste0(getwd(), "/Data for model/traps/", lyr = "traps_2B.shp"))

# mark-recapture data
load(paste0(getwd(), "/Data for model/MR/mr_data_2A.RData"))
load(paste0(getwd(), "/Data for model/MR/mr_data_2B.RData"))

# trap operation matrix
load(paste0(getwd(), "/Data for model/Trap op/trap_op_2A.RData"))
load(paste0(getwd(), "/Data for model/Trap op/trap_op_2B.RData"))

# ______________________________________________________________________________
# 4. Separate lists ----
# ______________________________________________________________________________
# 4a. Capture histories ----

# capture histories n x K
# we need to add NA columns for additional occasions (max = 10)

# ______________________________________________________________________________

# capture histories n x K
# we need to add NA columns for additional occasions (max = 10)
ch <- rbind(cbind(focal.mr.2A.2[[1]],
                  matrix(data = NA,
                         nrow = nrow(focal.mr.2A.2[[1]]),
                         ncol = 10 - ncol(focal.mr.2A.2[[1]]))),
            cbind(focal.mr.2B.2[[1]],
                  matrix(data = NA,
                         nrow = nrow(focal.mr.2B.2[[1]]),
                         ncol = 10 - ncol(focal.mr.2B.2[[1]]))))

# ______________________________________________________________________________
# 4b. Previous capture indicator ----

# same deal here

# ______________________________________________________________________________

prev.c <- rbind(cbind(focal.mr.2A.2[[2]],
                      matrix(data = NA,
                             nrow = nrow(focal.mr.2A.2[[2]]),
                             ncol = 10 - ncol(focal.mr.2A.2[[2]]))),
                cbind(focal.mr.2B.2[[2]],
                      matrix(data = NA,
                             nrow = nrow(focal.mr.2B.2[[2]]),
                             ncol = 10 - ncol(focal.mr.2B.2[[2]]))))

# ______________________________________________________________________________
# 4c. Individual indices ----
# ______________________________________________________________________________

# we;ll have to cross-reference this with the open CHs later
indivs <- 1:length(c(focal.mr.2A.2[[3]],
                     focal.mr.2B.2[[3]]))

# group index 
group <- c(rep(1, times = length(focal.mr.2A.2[[3]])),
           rep(2, times = length(focal.mr.2B.2[[3]])))


# sex indicator (0 is F)
sex <- ifelse(c(focal.mr.2A.2[[5]],
                focal.mr.2B.2[[5]]) == "F",
              0,
              1)

# ______________________________________________________________________________
# 5. Define required quantities ----

# quite a bit of this can be packed into a separate script later
# the S bounds, trap locations, etc. can all be pre-attributed
# especially since they only vary by g and not g x t

# ______________________________________________________________________________

# number of sessions
G = 2

# number of captured individuals (by session)
n.ind = c(nrow(focal.mr.2A.2[[1]]),
          nrow(focal.mr.2B.2[[1]]))

# number of traps (always 36)
J = 36

# number of occasions (by session!)
K = c(ncol(focal.mr.2A.2[[1]]),
      ncol(focal.mr.2B.2[[1]]))

# S limits - each [G x 2]
S.xlim <- matrix(data = NA, nrow = 2, ncol = 2)
S.ylim <- matrix(data = NA, nrow = 2, ncol = 2)

  # 2A
  S.xlim[1, ] = as.numeric(st_bbox(S.2A)[c(1, 3)])
  S.ylim[1, ] = as.numeric(st_bbox(S.2A)[c(2, 4)])
  
  # 2B
  S.xlim[2, ] = as.numeric(st_bbox(S.2B)[c(1, 3)])
  S.ylim[2, ] = as.numeric(st_bbox(S.2B)[c(2, 4)])

# trap coordinates - [J x 2 x G]
j.coords <- array(data = NA, dim = c(J, 2, G))

  # 2A
  j.coords[ , 1, 1] = st_coordinates(traps.2A)[ , 1]
  j.coords[ , 2, 1] = st_coordinates(traps.2A)[ , 2]
  
  # 2B
  j.coords[ , 1, 2] = st_coordinates(traps.2B)[ , 1]
  j.coords[ , 2, 2] = st_coordinates(traps.2B)[ , 2]

# ______________________________________________________________________________
# 6. Data augmentation ----
# ______________________________________________________________________________

n.aug <- vector(length = 2) 

n.aug[1] = n.ind[1] * 3
n.aug[2] = n.ind[2] * 3

M = sum(n.ind) + sum(n.aug)

# data augmentation variable
z <- c(rep(1, times = sum(n.ind)), rep(NA, times = sum(n.aug)))

# zeroes variable
zeroes <- c(rep(NA, times = sum(n.ind)), rep(0, times = sum(n.aug)))

# add capture histories (no need to do directly)
#ch.1 <- rbind(ch, matrix(J + 1,    # index for "no captures"
#                         nrow = M - sum(n.ind), 
#                         ncol = ncol(ch)))

# add to previous capture covariate
prev.c.1 <- rbind(prev.c, matrix(0,
                                 nrow = M - sum(n.ind), 
                                 ncol = ncol(prev.c)))

# the zeroes in the 9 and 10 columns are probably okay here since we're only 
# looping through K for each session (MIGHT HAVE TO CHANGE)

# group assignment
group.1 <- c(group,
             rep(1, times = n.aug[1]),
             rep(2, times = n.aug[2]))

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
  S.area = c(st_area(S.2A) * 0.0001,   # area in ha,
             st_area(S.2B) * 0.0001),
  
  # known indices
  K = K,
  group = group.1
  
)

# data
data.list <- list(
  
  # individual data
  ch = ch,
  prev.c = prev.c.1,
  sex = sex.1,
  
  # trap operation (WILL NEED TO MODIFY TO MAKE MATRICES CONFORM TO MAX(K))
  # okay to hard code for now
  trap.op = array(c(trap.op.2A, trap.op.2B), dim = c(J, 8, 2)),
  
  # state space
  S.xlim = S.xlim,
  S.ylim = S.ylim,
  
  # traps
  j.coords = j.coords,
  
  # data augmentation [M]
  z = z,
  zeroes = zeroes,
  
  # which groups to turn z on and off for N sum
  which.groups.1 = ifelse(group.1 == 1, 1, 0),
  which.groups.2 = ifelse(group.1 == 2, 1, 0)
  
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
  
  # data augmentation
  psi ~ dunif(0, 1)
  psi.sex ~ dunif(0, 1)
  
  # loop through all individuals [M]
  for (i in 1:M) {
    
    # data augmentation parameter
    z[i] ~ dbern(psi)
    
    # sex
    sex[i] ~ dbern(psi.sex)
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
  
  # density D
  D[1:G] <- N[1:G] / S.area[1:G]
  
})

# ______________________________________________________________________________
# 8b. Initial values ----

# 01-07-2026
# these should be random so each chain gets a different starting point

# ______________________________________________________________________________

inits <- list(
  
  # initial s - start at 0 ,0
  s = cbind(rep(0, times = M),
            rep(0, times = M)),
  
  sigma = c(25, 25),
  
  alpha0 = c(-2, -2),
  
  alpha2 = c(0, 0),
  
  psi = 0.5,
  
  psi.sex = 0.5,
  
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
  "D",
  "z"
  
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
  samplesAsCodaMCMC = TRUE
  
)
toc()

# ______________________________________________________________________________
# 9. Summarize and visualize ----
# ______________________________________________________________________________

# summary of main parameters
MCMCsummary(model.1.run,
            params = c("D"))

# traceplot
MCMCtrace(model.1.run,
          params = "psi",
          pdf = F)

# 01-07-2026
# well, I think I just about sped this up as much as I can. Now to run the whole thing!

