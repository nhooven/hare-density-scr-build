# PROJECT: Building hare SCR model
# SCRIPT: 01 - Data formatting
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 09 Dec 2025
# COMPLETED: 10 Dec 2025
# LAST MODIFIED: 10 Dec 2025
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
library(mefa4)

# ______________________________________________________________________________
# 3. Read in data ----
# ______________________________________________________________________________
# 3a. Open capture histories ----

# this includes all individuals and their mark-recapture IDs (MRID)

# ______________________________________________________________________________

open.ch <- read.csv("hare_open_CapHis.csv")

# ______________________________________________________________________________
# 3b. Mark-recapture data ----

# this is a J x K frame of MRIDs, that I will need to manipulate for SCR
# we'll extract both individual capture histories and a trap operation matrix from this

# ______________________________________________________________________________

mr.data <- read.csv("trap_data_2022.csv")

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
# 4b. Mark-recapture data ----

# here we will subset only 2A, pretty simple

# ______________________________________________________________________________

focal.mr <- mr.data %>%
  
  # filter site
  filter(Site == "2A") %>%
  
  # only keep columns with more than one value (i.e., true occasions)
  # note that this also drops the Site column
  dplyr::select(where(~ n_distinct(.) > 1))

# remove trap.id for 
focal.mr <- focal.mr[ , -1]

# ______________________________________________________________________________
# 4c. Check that lists of MRIDs from both match exactly ----

# we'll return a df including the following information:
# MRID
  # The MRIDs with discrepancies
# which.list
  # iMR.niIndivs = MRIDs in MR list but not in Indivs list
  # niMR.iIndivs = MRIDS not in MR list but in Indivs list
# where.else
  # which other unit MRIDs from MR were captured in during the same session (if applicable)

check_mrid <- function (
    
  indivs,
  mr,
  session
  
) {
  
  # generate list of MRIDs from entire trapping data
  # bind together into one vector
  all.outcomes <- as.vector(as.matrix(mr))
  
  # remove trap indicators
  all.outcomes.unique <- unique(all.outcomes)
  
  trap.mrid <- all.outcomes.unique[! all.outcomes.unique %in% c("O", "B", "X", "C", "E")]
  
  # compare to focal.indivs
  # which MRIDs are IN MR and NOTIN indivs?
  iMR.niIndivs <- trap.mrid[which(trap.mrid %notin% indivs$MRID)]
  
  # which MRIDs are NOTIN MR and IN indivs?
  niMR.iIndivs <- indivs$MRID[which(indivs$MRID %notin% trap.mrid)]
  
  # df to hold information
  check.df <- rbind(
    
    # find if those MRIDs were caught elsewhere during the same session
    data.frame(MRID = iMR.niIndivs,
               which.list = "iMR.niIndivs",
               where.else = open.ch[ , c(session)][which(open.ch$MRID %in% iMR.niIndivs)]),
    
    # list those in indivs but not in MR
    data.frame(MRID = ifelse(length(niMR.iIndivs) > 0,
                             niMR.iIndivs,
                             NA),
               which.list = "niMR.iIndivs",
               where.else = NA)
    
  )
  
  # return
  return(check.df)
  
}

# ______________________________________________________________________________

# apply function
(check.df <- check_mrid(focal.indivs, focal.mr, "MR22f"))

# ______________________________________________________________________________
# 4d. Change any MR entries from "other unit" individuals ----


# if we have individuals captured originally in (an)other unit(s),
# change their MR entry(ies) to "C"
# this is to account for the trap being partially closed during the night.
# We'll keep those individuals with their original unit and assume any
# their AC is somewhere between the units.

# an alternative is to merge the grids (for 2A/2B, for example)
# but this complicates any inference on treatment effects
# so for simplicity, since grid-switching (within year) is rare,
# we'll ignore it for the most part

# ______________________________________________________________________________

# function
change_entry_otherUnit <- function (
    
  check.df,
  mr
  
) {
  
  # define which MRIDs to switch
  MRIDs.to.switch <- check.df$MRID[check.df$which.list == "iMR.niIndivs"]

  # change to vector
  mr.vect <- as.vector(as.matrix(mr))
  
  # change to "C"
  mr.vect[which(mr.vect %in% MRIDs.to.switch)] <- "C"
  
  # convert back to df and return
  mr.2 <- as.data.frame(matrix(mr.vect, nrow = nrow(mr), ncol = ncol(mr)))
  
  colnames(mr.2) <- colnames(mr)
  
  return(mr.2)
  
}

# ______________________________________________________________________________

# use function
focal.mr.2 <- change_entry_otherUnit(check.df, focal.mr)

# ______________________________________________________________________________
# 4e. Reformat MR data ----

# now we need to cast this into the correct format for SCR
# a n x K matrix, where each entry is a trap number (from J)
# no capture by occasion K (NAs) will eventually be given the J + 1 index

# in this function, we'll also create the "trapped previously" binary matrix
# i.e., the trap response covariate!

# ______________________________________________________________________________

# function
castMR <- function(
  
  mr,
  indivs
  
) {
  
  # vector of individuals
  indivs.mrid <- indivs$MRID
  
  # set up matrix
  mr.mat <- matrix(data = NA,
                   nrow = length(indivs.mrid),
                   ncol = ncol(mr))
  
  # loop through MRIDs (i)
  for (i in 1:length(indivs.mrid)) {
    
    # loop through occasions (K)
    for (k in 1:ncol(mr)) {
      
      # find which trap (if any) the focal individual was found in [i, k]
      if (length(which(mr[, k] == indivs.mrid[i])) > 0) {
        
        mr.mat[i, k] <- which(mr[, k] == indivs.mrid[i])
        
        # if not, add the "not trapped" index (i.e., 37)
      } else {
        
        mr.mat[i, k] <- nrow(mr) + 1
        
      }
      
    }
    
  }
  
  # trap response covariate
  # binary - if previously captured, this will be zero
  prev.cap.mat <- matrix(data = NA,
                         nrow = length(indivs.mrid),
                         ncol = ncol(mr))
  
  # first occasion must be zero
  prev.cap.mat[ , 1] <- 0
  
  # loop through MRIDs (i)
  for (i in 1:length(indivs.mrid)) {
    
    # loop through occasions (2...K)
    for (k in 2:ncol(mr)) {
      
      # if was captured on ANY previous occasions, add a 1
      if (any(mr.mat[i, 1:(k - 1)] %in% 1:nrow(mr))) {
        
        prev.cap.mat[i, k] <- 1
        
        # if not, add a 0
      } else {
        
        prev.cap.mat[i, k] <- 0
        
      }
      
    }
    
  }
  
  # pack into list and return
  castMR.list <- list(
    
    mr.mat, 
    prev.cap.mat,
    1:length(indivs.mrid),        # individual index
    indivs.mrid
    
    )
  
  return(castMR.list)
  
}

# use function
focal.mr.3 <- castMR(focal.mr.2, focal.indivs)

# now these data are formatted correctly for SCR

# ______________________________________________________________________________
# 4f. Add sex covariate ----

# for now, the only covariate I'm reasonably sure we'll use is sex
# let's attribute from the focal.indivs frame

# ______________________________________________________________________________

focal.mr.3[[5]] <- focal.indivs$Sex

# ______________________________________________________________________________
# 5. Write to file ----

# we'll keep this as an .RData file for now to preserve the list structure (it's handy)

# ______________________________________________________________________________

save(focal.mr.3, file = paste0(getwd(), "/Data for model/mr_data.RData"))

# ______________________________________________________________________________
# 6. Trap operation matrix ----

# J x K matrix
# The values in here will be multiplied within the model to allow for or "zero out"
# trap-specific capture probabilities
# Cheekily, I will include 0.5 as a possible value to encode our uncertainty about
# trap availability given the trap was:
  # (1) found closed or 
  # (2) caught a bycatch species

# we won't know how long the trap was available for a hare capture, so we consider
# its availability to be a Bernoulli trial. Maybe it was open all night and closed
# or caught another critter in the morning after no hare was going to in. Maybe
# it closed or caught another critter shortly after we set it.

# In effect, the multiplication by 0.5 allows the full capture probability for
# half of all MCMC draws, and zeroes it out for the other half

# any "escaped" hares that did not receive an AnimalID will also induce a 0.5 here

# ______________________________________________________________________________

focal.mr.2.vect <- unlist(focal.mr.2)

focal.mr.2.vect[which(focal.mr.2.vect %notin% c("X", "C", "B", "E"))] <- 1
focal.mr.2.vect[which(focal.mr.2.vect %in% c("X"))] <- 0
focal.mr.2.vect[which(focal.mr.2.vect %in% c("C", "B", "E"))] <- 0.5

trap.op <- matrix(as.numeric(focal.mr.2.vect),
                  nrow = nrow(focal.mr.2),
                  ncol = ncol(focal.mr.2))

trap.op

# write to file
save(trap.op, file = paste0(getwd(), "/Data for model/trap_op.RData"))
