# PROJECT: Building hare SCR model
# SCRIPT: 01 - Data formatting
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 09 Dec 2025
# COMPLETED: 10 Dec 2025
# LAST MODIFIED: 08 Jan 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# In this script, I will input, clean, and format all the data I need to fit
# a closed SCR for all fall 2022 sessions (3B is NOT included)

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

# here we will subset only the individuals we need
# for the basic multi-session model, we'll pull 2A and 2B from 2022f

# function
subset_open_ch <- function (open.ch, session) {
  
  focal.indivs <- open.ch %>%
    
    filter(MR22f == session) %>%
    
    # keep AnimalID, MRID, and Sex
    dplyr::select(AnimalID.1,
                  MRID,
                  Sex)
  
  return(focal.indivs)
  
}

# ______________________________________________________________________________

# importantly, MRIDs need to match EXACTLY with their entries in the trapping data
# I'll implement this as a "check" function once everything is formatted correctly

focal.indivs.2A <- subset_open_ch(open.ch, session = "2A")
focal.indivs.2B <- subset_open_ch(open.ch, session = "2B")
focal.indivs.2C <- subset_open_ch(open.ch, session = "2C")
focal.indivs.3A <- subset_open_ch(open.ch, session = "3A")
focal.indivs.3C <- subset_open_ch(open.ch, session = "3C")

# ______________________________________________________________________________
# 4b. Mark-recapture data ----

# subset both units

# function
subset_mr <- function (mr.data, session) {
  
  focal.mr <- mr.data %>%
    
    # filter site
    filter(Site == session) %>%
    
    # only keep columns with more than one value (i.e., true occasions)
    # note that this also drops the Site column
    dplyr::select(where(~ n_distinct(.) > 1)) %>%
    
    # remove trap.id
    dplyr::select(-Trap)
  
  return(focal.mr)
  
}

# ______________________________________________________________________________

focal.mr.2A <- subset_mr(mr.data, "2A")
focal.mr.2B <- subset_mr(mr.data, "2B")
focal.mr.2C <- subset_mr(mr.data, "2C")
focal.mr.3A <- subset_mr(mr.data, "3A")
focal.mr.3C <- subset_mr(mr.data, "3C")

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
  
  # dfs to hold information
  # in MR and not in Indivs
  if (length(iMR.niIndivs) > 0) {
    
    iMR.niIndivs.df <- data.frame(
      
      MRID = iMR.niIndivs,
      which.list = "iMR.niIndivs",
      where.else = open.ch[ , c(session)][which(open.ch$MRID %in% iMR.niIndivs)]
      
      )
    
  } else {
    
    iMR.niIndivs.df <- 
      
      data.frame(
    
        MRID = NA,
        which.list = "iMR.niIndivs",
        where.else = NA
        
      )
    
  }
  
  # not in MR and in Indivs
  if (length(niMR.iIndivs) > 0) {
    
    niMR.iIndivs.df <- data.frame(
      
      MRID = niMR.iIndivs,
      which.list = "niMR.iIndivs",
      where.else = NA
      
    )
    
  } else {
    
    niMR.iIndivs.df <- data.frame(
      
      MRID = NA,
      which.list = "niMR.iIndivs",
      where.else = NA
      
    )
    
  }
  
  check.df <- rbind(iMR.niIndivs.df, niMR.iIndivs.df)
  
  # return
  return(check.df)
  
}

# ______________________________________________________________________________

# apply function
(check.df.2A <- check_mrid(focal.indivs.2A, focal.mr.2A, "MR22f"))
(check.df.2B <- check_mrid(focal.indivs.2B, focal.mr.2B, "MR22f"))
(check.df.2C <- check_mrid(focal.indivs.2C, focal.mr.2C, "MR22f"))
(check.df.3A <- check_mrid(focal.indivs.3A, focal.mr.3A, "MR22f"))
(check.df.3C <- check_mrid(focal.indivs.3C, focal.mr.3C, "MR22f"))

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
focal.mr.2A.1 <- change_entry_otherUnit(check.df.2A, focal.mr.2A)

# ______________________________________________________________________________
# 4e. Reformat MR data ----

# now we need to cast this into the correct format for SCR
# a n x K matrix, where each entry is a trap number (from J)
# no capture by occasion K (NAs) will eventually be given the J + 1 index

# in this function, we'll also create the "trapped previously" binary matrix
# i.e., the trap response covariate!

# this returns a convenient list:
  # [[1]]: capture histories formatted for categorical likelihood n x K
  # [[2]]: binary trap response matrix n x K
  # [[3]]: individual index n (should bind to site name so we can make this unique later)
  # [[4]]: MRID if we need it

# ______________________________________________________________________________

# function
castMR <- function(
  
  mr,
  indivs,
  site
  
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
    paste0(site, "_", 1:length(indivs.mrid)),        # individual index
    indivs.mrid
    
    )
  
  return(castMR.list)
  
}

# use function
focal.mr.2A.2 <- castMR(focal.mr.2A.1, focal.indivs.2A, "2A")
focal.mr.2B.2 <- castMR(focal.mr.2B, focal.indivs.2B, "2B")
focal.mr.2C.2 <- castMR(focal.mr.2C, focal.indivs.2C, "2C")
focal.mr.3A.2 <- castMR(focal.mr.3A, focal.indivs.3A, "3A")
focal.mr.3C.2 <- castMR(focal.mr.3C, focal.indivs.3C, "3C")

# now these data are formatted correctly for SCR

# ______________________________________________________________________________
# 4f. Add sex covariate ----

# for now, the only covariate I'm reasonably sure we'll use is sex
# let's attribute from the focal.indivs frame

# ______________________________________________________________________________

focal.mr.2A.2[[5]] <- focal.indivs.2A$Sex
focal.mr.2B.2[[5]] <- focal.indivs.2B$Sex
focal.mr.2C.2[[5]] <- focal.indivs.2C$Sex
focal.mr.3A.2[[5]] <- focal.indivs.3A$Sex
focal.mr.3C.2[[5]] <- focal.indivs.3C$Sex

# ______________________________________________________________________________
# 5. Write to file ----

# we'll keep these as .RData files for now to preserve the list structure (it's handy)

# ______________________________________________________________________________

save(focal.mr.2A.2, file = paste0(getwd(), "/Data for model/MR/mr_data_2A.RData"))
save(focal.mr.2B.2, file = paste0(getwd(), "/Data for model/MR/mr_data_2B.RData"))
save(focal.mr.2C.2, file = paste0(getwd(), "/Data for model/MR/mr_data_2C.RData"))
save(focal.mr.3A.2, file = paste0(getwd(), "/Data for model/MR/mr_data_3A.RData"))
save(focal.mr.3C.2, file = paste0(getwd(), "/Data for model/MR/mr_data_3C.RData"))

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

make_trap_op <- function (focal.mr) {
  
  focal.mr.vect <- unlist(focal.mr)

  focal.mr.vect[which(focal.mr.vect %notin% c("X", "C", "B", "E"))] <- 1
  focal.mr.vect[which(focal.mr.vect %in% c("X"))] <- 0
  focal.mr.vect[which(focal.mr.vect %in% c("C", "B", "E"))] <- 0.5
  
  trap.op <- matrix(as.numeric(focal.mr.vect),
                    nrow = nrow(focal.mr),
                    ncol = ncol(focal.mr))
  
  return(trap.op)
  
}

# use function
trap.op.2A <- make_trap_op(focal.mr.2A.1)
trap.op.2B <- make_trap_op(focal.mr.2B)
trap.op.2C <- make_trap_op(focal.mr.2C)
trap.op.3A <- make_trap_op(focal.mr.3A)
trap.op.3C <- make_trap_op(focal.mr.3C)

# write to file
save(trap.op.2A, file = paste0(getwd(), "/Data for model/Trap op/trap_op_2A.RData"))
save(trap.op.2B, file = paste0(getwd(), "/Data for model/Trap op/trap_op_2B.RData"))
save(trap.op.2C, file = paste0(getwd(), "/Data for model/Trap op/trap_op_2C.RData"))
save(trap.op.3A, file = paste0(getwd(), "/Data for model/Trap op/trap_op_3A.RData"))
save(trap.op.3C, file = paste0(getwd(), "/Data for model/Trap op/trap_op_3C.RData"))

