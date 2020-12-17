rm(list=ls())
library(deSolve)
library(tidyverse)
library(lubridate)
library(reshape2)
library(parallel)
library(scales)
library(nonpar)
set.seed(1211)

source("sd_funcs.R") # this version has sd=0.15 on R0

if(file.exists('Figure2.RData')){load('Figure2.RData')}else{
  
  # initialize DE model
  N = 5.1e6  #population of BC
  pars = list(N = N,
              D = 5,
              R0 = 2.57, # from pre-distancing fit
              k1 = 1/5,
              k2 = 1,
              q = 0.05,
              r = 1,
              ur = 0.4,
              f = 1.0) # initially no distancing
  
  fsi = with(pars,
             r/(r+ur))
  nsi = 1 - fsi
  i0 = 40
  # initial state of the ode model
  state = c(S = nsi*(N-i0),
            E1 = 0.4*nsi*i0,
            E2 = 0.1*nsi*i0,
            I = 0.5*nsi*i0,
            Q = 0,
            R = 0,
            Sd = fsi*(N-i0),
            E1d = 0.4*fsi*i0,
            E2d = 0.1*fsi*i0,
            Id = 0.5*fsi*i0,
            Qd = 0,
            Rd = 0)
  
  # Number of R0s to sample
  nReps = 100
  
  time_list = list()
  longtimes = seq(from = -20,
                  to = 120,
                  by = 0.1) # how long to simulate for
  
  
  # general function for social distancing timing
  timing = function( t, f1, f2, tstar ){
    ifelse( t<tstar, f1, f2)
  }
  
  
  # Start on March 1st, day 0
  start_date = dmy("01-03-2020")
  
  # social distancing started March 18th (day 17)
  sd_date = dmy("18-03-2020")
  
  # number of days between start and SD implemented
  t_sd = as.numeric(sd_date)-as.numeric(start_date)
  
  
  # No distancing : f stays constant at 1.0
  time_list[[1]] = function(t, pars) {timing(t, pars$f, pars$f, t_sd)} # always 1.0
  
  # Strong dist : f -> 0.3 at sd_date
  time_list[[2]] = function(t, pars) {timing(t, pars$f, 0.3, t_sd)}
  
  # Medium dist : f -> 0.4 at sd_date
  time_list[[3]] = function(t, pars) {timing(t, pars$f, 0.4, t_sd)}
  
  # Weak dist : f => 0.7 at sd_date
  time_list[[4]] = function(t, pars) {timing(t, pars$f, 0.7, t_sd)}
  
  
  # generate trajectories
  # uncertainty arises from our uncertainty in our estimate of R0
  tt_no_noise = lapply(time_list, function(x) multisolve(params=pars,
                                                         state,
                                                         longtimes,
                                                         nReps = nReps,
                                                         timefunc=x))
  names(tt_no_noise) = c(1.0, 0.3, 0.4, 0.7)
  
  # Calculate the first time the curves diverge
  firstdiv = firstdiverge2(tt_no_noise, increasing=F)
  
  save.image('Figure2.RData')
}

plottimes = seq(from = -20,
                to = 100,
                by = 0.1)

rbind(cbind( getCasesbyDay2(tt_no_noise[[1]],plottimes,startDate=start_date), f2="baseline"),
      cbind( getCasesbyDay2(tt_no_noise[[3]],plottimes,startDate=start_date), f2="0.4"),
      cbind( getCasesbyDay2(tt_no_noise[[4]],plottimes,startDate=start_date), f2="0.7")) %>%
  ggplot(aes(x=dates, y=median))+
  geom_line(aes(colour=f2))+
  geom_ribbon(aes(x=dates, ymin=lower25, ymax=upper75, fill=f2), alpha=0.3)+
  scale_fill_manual(values=c("0.4"="#7ec5a3", "0.7"="#f7bb82", "baseline"="#cdcdcd"))+
  scale_colour_manual(values=c("0.4"="#7ec5a3", "0.7"="#f7bb82", "baseline"="#cdcdcd"))+
  geom_vline(xintercept=start_date+firstdiv$results$time[2:3], colour=c("#7ec5a3", "#f7bb82"), size=1, linetype="dashed")+
  geom_vline(xintercept=dmy("18-03-2020"), colour="#cdcdcd", size=1)+
  xlim(dmy("01-03-2020"), dmy("15-04-2020"))+
  coord_cartesian(ylim = c(0,6000)) +
  labs(x="Date", y="Active cases", fill=expression(f[1]), colour=expression(f[1]))+
  theme_minimal()

p = last_plot()
p
ggsave('Figure2a.pdf', width=6, height=2, units = 'in', dpi = 900)
