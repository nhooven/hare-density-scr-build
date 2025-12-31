# PROJECT: Building hare SCR model
# SCRIPT: 05 - Goodness-of-fit tests
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 31 Dec 2025
# COMPLETED: 31 Dec 2025
# LAST MODIFIED: 31 Dec 2025
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# Here we'll calculate each fit statistic post-hoc
# this is mostly because I'm unsure how to do this
# maybe eventually we can calculate it within the model

# ______________________________________________________________________________
# 2. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(coda)
library(tictoc)
library(MCMCvis)
library(cowplot)

# ______________________________________________________________________________
# 3. Read in data ----
# ______________________________________________________________________________

load(paste0(getwd(), "/Model draws/model_1.RData"))
load(paste0(getwd(), "/Model data lists/const_1.RData"))
load(paste0(getwd(), "/Model data lists/data_1.RData"))

# ______________________________________________________________________________
# 4. Goodness-of-fit tests ----

# importantly, we should monitor all parameters from the model to do this
# we could save p[i, j, k] but eventually this would lead to a prohibitively
# large output file. We can calculate these probs easily as long as we save
# the relevant parameters

# ______________________________________________________________________________
# 4a. Define function ----

# this returns three plots showing the discrepancies

# ______________________________________________________________________________

SCR_GoF <- function (
    
  draws,
  iter = 10,
  data.list,
  const.list = constant.list
  
) {
  
  # model draws
  draws.1 = draws[sample(1:nrow(draws), size = iter), ]
  
  # loop through iterations (possibly a sample)
  FT.indiv.all <- data.frame()
  FT.trap.all <- data.frame()
  FT.indiv.trap.all <- data.frame()
  
  for (n in 1:nrow(draws.1)) {
    
    # calculate alpha1
    alpha1 <- vector(length = 2)
    
    alpha1[1] = -1 / draws.1[n, "sigma[1]"]
    alpha1[2] = -1 / draws.1[n, "sigma[2]"]
    
    # calculate trap-specific probabilities using parameter draws
    # loop through individuals
    s <- matrix(data = NA, nrow = const.list$M, ncol = 2)
    sex2 <- vector(length = const.list$M)
    
    # arrays to hold probabilities
    lp <- array(data = NA, dim = c(const.list$M, const.list$K, const.list$J))
    p <- array(data = NA, dim = c(const.list$M, const.list$K, const.list$J + 1))
    
    for (i in 1:const.list$M) {
      
      # activity center
      s[i, 1] = draws.1[n, paste0("s[", i, ", 1]")]
      s[i, 2] = draws.1[n, paste0("s[", i, ", 2]")]
      
      # sex indicator
      sex2[i] = draws.1[n, paste0("sex[", i, "]")] + 1
      
    # calculate distances
    # loop through traps
    d <- matrix(data = NA, nrow = const.list$M, ncol = const.list$J)
    
    for (j in 1:const.list$J) {
        
      # distance from each AC s to each trap j
      d[i, j] <- sqrt((s[i, 1] - data.list$x.j[j])^2 + (s[i, 2] - data.list$y.j[j])^2)
      
    }
    
    # calculate linear predictor FIRST
    for (k in 1:const.list$K) {
      
      for (j in 1:const.list$J) {
        
        # linear predictor (with trap operation matrix)
        # importantly, we won't multiply by z because we want the simulation to be naive to whether the individual is included
        lp[i, k, j] <- exp(draws.1[n, paste0("alpha0[", sex2[i], "]")] + 
                             
                             alpha1[sex2[i]] * d[i, j] + 
                             
                             draws.1[n, paste0("alpha2[", sex2[i], "]")] * data.list$prev.c[i, k]) * data.list$trap.op[j, k]
        
      }
        
    }
    
    # calculate probabilities
    for (k in 1:const.list$K) {
      
      for (j in 1:const.list$J) {
        
        p[i, k, j] <- lp[i, k, j] / (1 + sum(lp[i, k, 1:const.list$J]))  # sum over all traps
        
        # probability of not being captured as the complement of all trap-specific probs
        p[i, k, const.list$J + 1] <- 1 - sum(p[i, k, 1:const.list$J])
    
      }
      
    }
    
    }
    
    # simulate capture histories
    # these will be categorical- then we need to expand into a binary array
    ch.new <- matrix(data = NA, nrow = const.list$M, ncol = const.list$K)
    ch.new.array <- array(data = NA, dim = c(const.list$M,
                                             const.list$K,
                                             const.list$J))
    
    for (i in 1:const.list$M) {
      
      for (k in 1:const.list$K) {
        
        ch.new[i, k] <- rcat(prob = p[i, k, ])
        
        # expand into binary array
        # which cap
        binary.vect <- rep(0, times = const.list$J)
        
        if (ch.new[i, k] == const.list$J + 1) {
          
          ch.new.array[i, k, ] <- binary.vect
          
        } else {
          
          binary.vect[ch.new[i, k]] <- 1
          
          ch.new.array[i, k, ] <- binary.vect
          
        }
        
      }
      
    }
    
    # expand true capture histories
    ch.array <- array(data = NA, dim = c(const.list$M,
                                         const.list$K,
                                         const.list$J))
    
    for (i in 1:nrow(data.list$ch)) {
      
      for (k in 1:const.list$K) {
      
      # which cap
      binary.vect <- rep(0, times = const.list$J)
      
      if (data.list$ch[i, k] == const.list$J + 1) {
        
        ch.array[i, k, ] <- binary.vect
        
      } else {
        
        binary.vect[data.list$ch[i, k]] <- 1
        
        ch.array[i, k, ] <- binary.vect
        
      }
      
    }
        
  }
    
    # aggregate capture histories for test statistics
    
    # FT-indiv
    # these are individual encounter frequencies
    FT.indiv <- data.frame(
        
        # simulated counts
        count.sim = apply(ch.new.array, 1, sum),
        
        # real counts
        count.real = apply(ch.array, 1, sum),
        
        # expected counts (must remove the "no capture" bin)
        E.count = apply(p[ , , -(const.list$J + 1)], 1, sum)
        
        ) %>%
        
        mutate(
          
          err.sim = (sqrt(count.sim) - sqrt(E.count))^2,
          err.real = (sqrt(count.real) - sqrt(E.count))^2
          
        ) %>%
        
        summarize(
          
          err.sim.sum = sum(err.sim),
          err.real.sum = sum(err.real)
          
        )
      
      # bind in
      FT.indiv.all <- rbind(FT.indiv.all, FT.indiv)
    
    # FT-trap
    # these are trap-specific counts
    FT.trap <- data.frame(
        
        # simulated counts
        count.sim = apply(ch.new.array, 3, sum),
        
        # real counts
        count.real = apply(ch.array, 3, sum),
        
        # expected counts (must remove the "no capture" bin)
        E.count = apply(p[ , , -(const.list$J + 1)], 3, sum)
        
      ) %>%
        
        mutate(
          
          err.sim = (sqrt(count.sim) - sqrt(E.count))^2,
          err.real = (sqrt(count.real) - sqrt(E.count))^2
          
        ) %>%
        
        summarize(
          
          err.sim.sum = sum(err.sim),
          err.real.sum = sum(err.real)
          
        )
      
      # bind in
      FT.trap.all <- rbind(FT.trap.all, FT.trap)
    
    # FT-indiv-trap
    # these are individual-trap counts
    # should have low power according to Royle et al. 
    FT.indiv.trap <- list(
        
        # simulated counts
        count.sim = apply(ch.new.array, c(1, 3), sum),
        
        # real counts
        count.real = apply(ch.array, c(1, 3), sum),
        
        # expected counts (must remove the "no capture" bin)
        E.count = apply(p[ , , -(const.list$J + 1)], c(1, 3), sum)
        
      )
      
      # calculate by indiv-trap
      err.sim <- matrix(data = NA, nrow = const.list$M, ncol = const.list$J)
      err.real <- matrix(data = NA, nrow = const.list$M, ncol = const.list$J)
      
      for (i in 1:const.list$M) {
        
        for (j in 1:const.list$J) {
          
          err.sim[i, j] <- (sqrt(FT.indiv.trap$count.sim[i, j]) - sqrt(FT.indiv.trap$E.count[i, j]))^2
          
          err.real[i, j] <- (sqrt(FT.indiv.trap$count.real[i, j]) - sqrt(FT.indiv.trap$E.count[i, j]))^2
          
        }
        
      }
      
      FT.indiv.trap.df <- data.frame(err.sim.sum = sum(err.sim),
                                     err.real.sum = sum(err.real))
      
      # bind in
      FT.indiv.trap.all <- rbind(FT.indiv.trap.all, FT.indiv.trap.df)
    
  }
  
  # return plots and summarize
  # Bayesian p-values
  FT.indiv.test <- mean(FT.indiv.all$err.sim.sum > FT.indiv.all$err.real.sum)
  FT.trap.test <- mean(FT.trap.all$err.sim.sum > FT.trap.all$err.real.sum)
  FT.indiv.trap.test <- mean(FT.indiv.trap.all$err.sim.sum > FT.indiv.trap.all$err.real.sum)
  
  # plots
  FT.indiv.plot <- ggplot(data = FT.indiv.all) +
    
    theme_bw() +
    
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dashed") +
    
    geom_point(aes(x = err.real.sum,
                   y = err.sim.sum),
               shape = 21,
               color = "darkblue") +
    
    coord_cartesian(xlim = c(min(min(FT.indiv.all$err.sim.sum), min(FT.indiv.all$err.real.sum)),
                             max(max(FT.indiv.all$err.sim.sum), max(FT.indiv.all$err.real.sum))),
                    ylim = c(min(min(FT.indiv.all$err.sim.sum), min(FT.indiv.all$err.real.sum)),
                             max(max(FT.indiv.all$err.sim.sum), max(FT.indiv.all$err.real.sum)))) +
      
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black")) +
      
      xlab("Observed discrepancy") +
      
      ylab("Simulated discrepancy") +
      
      ggtitle(paste0("FT-individual - ", round(FT.indiv.test, digits = 2)))
  
  # FT-traps
  FT.trap.plot <- ggplot(data = FT.trap.all) +
    
    theme_bw() +
    
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dashed") +
    
    geom_point(aes(x = err.real.sum,
                   y = err.sim.sum),
               shape = 21,
               color = "darkblue") +
    
    coord_cartesian(xlim = c(min(min(FT.trap.all$err.sim.sum), min(FT.trap.all$err.real.sum)),
                             max(max(FT.trap.all$err.sim.sum), max(FT.trap.all$err.real.sum))),
                    ylim = c(min(min(FT.trap.all$err.sim.sum), min(FT.trap.all$err.real.sum)),
                             max(max(FT.trap.all$err.sim.sum), max(FT.trap.all$err.real.sum)))) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black")) +
    
    xlab("Observed discrepancy") +
    
    ylab("Simulated discrepancy") +
    
    ggtitle(paste0("FT-trap - ", round(FT.trap.test, digits = 2)))
  
  # FT-individual-traps
  FT.indiv.trap.plot <- ggplot(data = FT.indiv.trap.all) +
    
    theme_bw() +
    
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dashed") +
    
    geom_point(aes(x = err.real.sum,
                   y = err.sim.sum),
               shape = 21,
               color = "darkblue") +
    
    coord_cartesian(xlim = c(min(min(FT.indiv.trap.all$err.sim.sum), min(FT.indiv.trap.all$err.real.sum)),
                             max(max(FT.indiv.trap.all$err.sim.sum), max(FT.indiv.trap.all$err.real.sum))),
                    ylim = c(min(min(FT.indiv.trap.all$err.sim.sum), min(FT.indiv.trap.all$err.real.sum)),
                             max(max(FT.indiv.trap.all$err.sim.sum), max(FT.indiv.trap.all$err.real.sum)))) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black")) +
    
    xlab("Observed discrepancy") +
    
    ylab("Simulated discrepancy") +
    
    ggtitle(paste0("FT-indiv.trap - ", round(FT.indiv.trap.test, digits = 2)))

  # add together  
  all.plots <- plot_grid(FT.indiv.plot, FT.trap.plot, FT.indiv.trap.plot,
                         nrow = 1)

  return(all.plots)
  
}

# ______________________________________________________________________________
# 4b. Use function ----
# ______________________________________________________________________________

SCR_GoF(draws = model.1.run, iter = 1000, data.list = data.list)
