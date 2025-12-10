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
  occ            # total, all, or specific occasions
  
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
  # for now, no need to spatialize any of this
  # maybe if I want to include the unit boundary later or something I can
  
  # bounding box for traps
  traps.bbox <- st_bbox(traps.sf)
  
  # total - all in one map
  if (occ == "total") {
    
    # add in coords
    total.forPlot <- trap.cap.freq.all %>%
      
      mutate(x = st_coordinates(traps.sf)[ , 1],
             y = st_coordinates(traps.sf)[ , 2])
    
    # plot
    out.plot <- ggplot() +
      
      theme_bw() +
      
      geom_point(data = total.forPlot,
                 aes(x = x,
                     y = y,
                     fill = freq),
                 size = 4,
                 shape = 21,
                 stroke = 1) +
      
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            legend.title = element_blank()) +
      
      scale_fill_gradient(low = "white", high = "darkgreen") +
      
      # expand the plotting window a bit
      coord_cartesian(xlim = c(traps.bbox[1] - 40,
                               traps.bbox[3] + 40),
                      ylim = c(traps.bbox[2] - 40,
                               traps.bbox[4] + 40))
    
  }
  
  # add in coords
  all.forPlot <- trap.cap.freq.occ %>%
    
    mutate(x = st_coordinates(traps.sf)[ , 1],
           y = st_coordinates(traps.sf)[ , 2]) %>%
    
    # pivot for plotting
    pivot_longer(cols = 2:(ncol(mr) + 1)) %>%
    
    # factors
    mutate(name = factor(name,
                         labels = paste0("k == ", unique(name))))
  
  # "all" - facetted by occasion
  if (occ == "all") {
    
    # facetted plot
    # plot
    out.plot <- ggplot(all.forPlot) +
      
      theme_bw() +
      
      facet_wrap(~ name) +
      
      geom_point(aes(x = x,
                     y = y,
                     fill = as.factor(value)),
                 size = 2.5,
                 shape = 21,
                 stroke = 0.25) +
      
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            strip.background = element_rect(fill = "white"),
            strip.text = element_text(hjust = 0)) +
      
      scale_fill_manual(values = c("white", "darkgreen")) +
      
      # expand the plotting window a bit
      coord_cartesian(xlim = c(traps.bbox[1] - 40,
                               traps.bbox[3] + 40),
                      ylim = c(traps.bbox[2] - 40,
                               traps.bbox[4] + 40))
    
  }
  
  if (occ %in% 1:ncol(mr)) {
    
    # plot
    out.plot <- ggplot(all.forPlot %>% filter(name == paste0("k == ", occ))) +
      
      theme_bw() +
      
      geom_point(aes(x = x,
                     y = y,
                     fill = as.factor(value)),
                 size = 4,
                 shape = 21,
                 stroke = 1) +
      
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none") +
      
      scale_fill_manual(values = c("white", "darkgreen")) +
      
      # expand the plotting window a bit
      coord_cartesian(xlim = c(traps.bbox[1] - 40,
                               traps.bbox[3] + 40),
                      ylim = c(traps.bbox[2] - 40,
                               traps.bbox[4] + 40)) +
      
      ggtitle(paste0("k == ", occ))
    
  }
  
  return(out.plot)
  
}

total_caps_byTrap(focal.mr.3[[1]], traps.sf, occ = "total")
total_caps_byTrap(focal.mr.3[[1]], traps.sf, occ = "all")
total_caps_byTrap(focal.mr.3[[1]], traps.sf, occ = 2)

# ______________________________________________________________________________
# 5. Individual spatial capture histories ----
# ______________________________________________________________________________

# function
indiv_ch <- function (
    
  mr,
  traps.sf,
  indiv                # eventually it would be nice to allow for index by MRID
  
) {
  
  # subset mr
  focal.mr <- mr[indiv, ]
  
  # tabulate
  indiv.tab <- as.data.frame(table(focal.mr))
  
  # convert to integer
  indiv.tab$focal.mr <- as.integer(as.character(indiv.tab$focal.mr))
  
  # capture history
  focal.ch <- rbind(
    
    indiv.tab,
    
    data.frame(focal.mr = which(1:nrow(traps.sf) %notin% indiv.tab$focal.mr),
               Freq = 0)
    
  ) %>%
    
    # rename
    rename(trap.id = focal.mr) %>%
    
    # arrange by trap
    arrange(trap.id) %>%
    
    # remove "no captures"
    filter(trap.id != nrow(traps.sf) + 1) %>%
    
    mutate(
      
      # add coordinates
      x = st_coordinates(traps.sf)[ , 1],
      y = st_coordinates(traps.sf)[ , 2],
      
      # and add MRID
      MRID = focal.mr.3[[4]][indiv],
      
      # and sex
      sex = focal.mr.3[[5]][indiv])
  
  # plot
  # bounding box for traps
  traps.bbox <- st_bbox(traps.sf)
  
  # plot
  out.plot <- ggplot() +
    
    theme_bw() +
    
    geom_point(data = focal.ch,
               aes(x = x,
                   y = y,
                   fill = as.factor(Freq)),
               size = 4,
               shape = 21,
               stroke = 1) +
    
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.title = element_blank()) +
    
    scale_fill_viridis_d(option = "mako",
                         direction = -1) +
    
    # expand the plotting window a bit
    coord_cartesian(xlim = c(traps.bbox[1] - 40,
                             traps.bbox[3] + 40),
                    ylim = c(traps.bbox[2] - 40,
                             traps.bbox[4] + 40)) +
    
    ggtitle(paste0(focal.ch$MRID[1], " - ", focal.ch$sex[1]))
    
    
  return(out.plot)
  
}

# try it out
indiv_ch(focal.mr.3[[1]], traps.sf, 29)

# ______________________________________________________________________________
# 6. All spatial capture histories ----

# we could include both jittered capture locations and paths (a la "secr" package)

# ______________________________________________________________________________

# function
all_ch <- function (
    
  mr,
  traps.sf
  
) {
  
  # bounding box for traps
  traps.bbox <- st_bbox(traps.sf)
  
  # trap coords
  traps.df <- data.frame(x = st_coordinates(traps.sf)[ , 1],
                         y = st_coordinates(traps.sf)[ , 2])
  
  # loop through individuals
  all.ch <- data.frame()
  
  for (i in 1:nrow(mr)) {
    
  # subset mr
  focal.mr <- mr[i, ]
    
  # tabulate
  indiv.tab <- as.data.frame(table(focal.mr))
  
  # convert to integer
  indiv.tab$focal.mr <- as.integer(as.character(indiv.tab$focal.mr))
  
  # capture history
  focal.ch <- rbind(
    
    indiv.tab,
    
    data.frame(focal.mr = which(1:nrow(traps.sf) %notin% indiv.tab$focal.mr),
               Freq = 0)
    
  ) %>%
    
    # rename
    rename(trap.id = focal.mr) %>%
    
    # arrange by trap
    arrange(trap.id) %>%
    
    # remove "no captures"
    filter(trap.id != nrow(traps.sf) + 1) %>%
    
    mutate(
      
      # add coordinates
      x = st_coordinates(traps.sf)[ , 1],
      y = st_coordinates(traps.sf)[ , 2],
      
      # and add MRID
      MRID = focal.mr.3[[4]][i],
      
      # and sex
      sex = focal.mr.3[[5]][i]) %>%
    
    # replicate rows for each capture 
    slice(rep(1:n(), times = c(Freq))) %>%
    
    # and filter out traps without captures
    filter(Freq != 0) %>%
      
    # and drop Freq
    dplyr::select(-Freq)
  
  # and bind in
  all.ch <- rbind(all.ch, focal.ch)
  
  }
  
  
  # plot
  out.plot <- ggplot() +
    
    theme_bw() +
    
    # trap points as pluses
    geom_point(data = traps.df,
               aes(x = x,
                   y = y), 
               shape = 3,
               size = 2) +
    
    # paths show spatial recaptures
    geom_path(data = all.ch,
              aes(x = x,
                  y = y,
                  group = MRID,
                  color = MRID,
                  linetype = sex),
              linewidth = 1.5) +
    
    # jittered points so you can actually see the relative intensity of recaps
    geom_jitter(data = all.ch,
                aes(x = x,
                    y = y,
                    fill = MRID),
                size = 2,
                shape = 21,
                width = 7,
                height = 7) +
    
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "none") +
    
    # expand the plotting window a bit
    coord_cartesian(xlim = c(traps.bbox[1] - 40,
                             traps.bbox[3] + 40),
                    ylim = c(traps.bbox[2] - 40,
                             traps.bbox[4] + 40))
  
  
  return(out.plot)
  
}

# try it (F = solid, M = dashed)
all_ch(focal.mr.3[[1]], traps.sf)



# 12-10-2025
# I'll keep this here for now. I think I'm ready to fit some models tomorrow!