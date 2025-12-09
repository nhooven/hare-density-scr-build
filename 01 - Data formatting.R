# PROJECT: Building hare SCR model
# SCRIPT: 01 - Data formatting
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 09 Dec 2025
# COMPLETED: 
# LAST MODIFIED: 09 Dec 2025
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# In this script, I will input, clean, and format all the data I need to fit
# a closed SCR for one session (site 2A in fall 2022)

# I'll build models from the ground up, and eventually implement multiple sessions
# and an open population structure to estimate apparent survival and recruitment

# A few relevant variable conventions that I will adhere to throughout this process:

# n - total individuals captured
# J - total traps
# K - total occasions per session

# "focal" will refer to anything related to the one session of interest

# ______________________________________________________________________________
# 2. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 3. Read in data ----
# ______________________________________________________________________________
# 3a. Open capture histories ----

# this includes all individuals and their mark-recapture IDs (MRID)

# ______________________________________________________________________________

open.ch <- read.csv("hare_open_CapHis.csv")

# ______________________________________________________________________________
# 3b. Trapping data ----

# this is a J x K frame of MRIDs, that I will need to manipulate for SCR
# we'll extract both individual capture histories and a trap operation matrix from this

# ______________________________________________________________________________

trapping.data <- read.csv("trap_data_2022.csv")

# ______________________________________________________________________________
# 4. Clean data ----
# ______________________________________________________________________________
# 4a. "Open" capture histories ----

# here we will subset only the individuals we need (indivs caught during 2A 2022f)

# ______________________________________________________________________________

focal.indivs <- open.ch %>%
  
  # filter 2A during 2022f only
  filter(MR22f == "2A") %>%
  
  # keep AnimalID, MRID, and Sex
  dplyr::select(AnimalID.1,
                MRID,
                Sex)

# importantly, MRIDs need to match EXACTLY with their entries in the trapping data
# I'll implement this as a "check" function once everything is formatted correctly

# ______________________________________________________________________________
# 4b. Trapping data ----

# here we will subset only 2A, pretty simple

# ______________________________________________________________________________

focal.trapping <- trapping.data %>%
  
  # filter site
  filter(Site == "2A") %>%
  
  # only keep columns with more than one value (i.e., true occasions)
  # note that this also drops the Site column
  dplyr::select(where(~ n_distinct(.) > 1))

# ______________________________________________________________________________
# 4c. Check that lists of MRIDs from both match exactly ----
# ______________________________________________________________________________

check_mrid <- function (
    
  indivs,
  trapping
  
) {
  
  # generate list of MRIDs from entire trapping data
  # bind together into one vector
  all.outcomes <- as.vector(as.matrix(trapping[ , c(2:ncol(trapping))]))
  
  # remove trap indicators
  all.outcomes.unique <- unique(all.outcomes)
  
  trap.mrid <- all.outcomes.unique[! all.outcomes.unique %in% c("O", "B", "X", "C", "E")]
  
  # compare to 
  
}




# reformat trapping data
# now we need to cast this into the correct format for SCR
# a n x K matrix, where each entry is a trap number (from J)
# no capture by occasion K (NAs) will eventually be given the J + 1 index

focal.cap.his <- focal.trapping %>%
  
  pivot_longer(cols = D1:last_col(),
               names_to = "K")



