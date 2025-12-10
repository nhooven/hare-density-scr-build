# PROJECT: Building hare SCR model
# SCRIPT: 03 - Data exploration
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 10 Dec 2025
# COMPLETED: 
# LAST MODIFIED: 10 Dec 2025
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# It will be helpful to examine spatial capture histories, as well as total
# spatial captures, as well as summarizing various information on recaptures, etc.

# ______________________________________________________________________________
# 2. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(mefa4)

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
# 4. Total captures by trap ----

# We'll allow this to vary by occasion, or look at the total

# ______________________________________________________________________________

# function
total_caps_byTrap <- function (
    
  mr,
  traps.sf,
  occ = "all"
  
) {
  
  # compute totals by occasion
  total.list <- apply(mr, 2, table)
  
  all.occ <- data.frame()
  
  for (i in 1:length(total.list)) {
    
    focal.occ <- as.data.frame(total.list[[i]])
    
    # coerce to integer
    focal.occ$Var1 <- as.integer(as.character(focal.occ$Var1))
    
    # add traps that did not catch anything
    focal.occ.1 <- rbind(
      
      focal.occ,
      
      data.frame(Var1 = which(1:nrow(traps.sf) %notin% focal.occ$Var1),
                 Freq = 0)
      
    ) %>%
      
      # rename
      rename(trap.id = Var1) %>%
      
      # arrange by trap
      arrange(trap.id) %>%
      
      # remove "no captures"
      filter(trap.id != nrow(traps.sf) + 1) %>%
      
      # add occasion variable
      mutate(occ = i)
    
    # bind in
    all.occ <- rbind(all.occ, focal.occ.1)
    
  }
  
  # pivot - by occasion
  trap.cap.freq.occ <- all.occ %>%
    
    pivot_wider(names_from = occ,
                values_from = Freq)
  
  # summarize in total
  trap.cap.freq.all <- data.frame(trap.id = trap.cap.freq.occ$trap.id,
                                  freq = apply(trap.cap.freq.occ[ , c(2:ncol(trap.cap.freq.occ))], 1, sum))
  
  # plot
  if (occ == "all") {
    
    
    
  }
  
}