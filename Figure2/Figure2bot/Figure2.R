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
  # read in data
  bcdata = read.csv("../../bc_casecounts2204.csv", header=TRUE)

  bcdata[,1] = dmy(bcdata[,1])
  bcdata = bcdata[41:53,]
  bcdata$Diffs = c(bcdata$Cases[2] - bcdata$Cases[1],
                   diff(bcdata$Cases))
  bcdata$day = 1:nrow(bcdata)

  # initialize DE model
  N = 5.1e6  #population of BC
  pars = list(N = N,
              D = 5,
              R0 = 2.57,
              k1 = 1/5,
              k2 = 1,
              q = 0.05,
              r = 1,
              ur = 0.4,
              f = 0.36)

  fsi = with(pars,
             r/(r+ur))
  nsi = 1 - fsi
  i0 = 470

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

  # simulation
  nReps = 100

  time_list = list()
  longtimes = seq(from = -20,
                  to = 120,
                  by = 0.1) # how long to simulate for


  # general function for social distancing timing
  timing = function( t, f1, f2, tstar ){
    ifelse( t<tstar, f1, f2)
  }


  # Begin on April 10th
  start_date = dmy("10-04-2020")

  # social distancing relaxed starting May 17th
  relax_date = dmy("17-05-2020")

  # number of days between social distancing implemented and relaxed
  t_relax = as.numeric(relax_date)-as.numeric(start_date)


  # No relaxation of distancing : f stays constant at 0.36
  time_list[[1]] = function(t, pars) {timing(t, pars$f, pars$f, t_relax)} # always 0.36

  # Weak relaxation : f -> 0.5 at relax_date
  time_list[[2]] = function(t, pars) {timing(t, pars$f, 0.5, t_relax)}

  # Medium relaxation : f -> 0.65 at relax_date
  time_list[[3]] = function(t, pars) {timing(t, pars$f, 0.65, t_relax)}

  # Strong relaxation : f => 0.9 at relax_date
  time_list[[4]] = function(t, pars) {timing(t, pars$f, 0.9, t_relax)}

  # generate trajectories with no measurement noise
  # uncertainty arises from our uncertainty in our estimate of R0
  tt_no_noise = lapply(time_list, function(x) multisolve(params=pars,
                                                         state,
                                                         longtimes,
                                                         nReps = nReps,
                                                         timefunc=x))
  names(tt_no_noise) = c(0.36, 0.5, 0.65, 0.9)

  firstdiv = firstdiverge2(tt_no_noise, increasing=T)

  save.image('Figure2.RData')
}

plottimes = seq(from = -20,
                to = 100,
                by = 0.1)

rbind(cbind( getCasesbyDay2(tt_no_noise[[1]],plottimes,startDate=start_date), f2="baseline"),
      cbind( getCasesbyDay2(tt_no_noise[[2]],plottimes,startDate=start_date), f2="0.50"),
      cbind( getCasesbyDay2(tt_no_noise[[3]],plottimes,startDate=start_date), f2="0.65"),
      cbind( getCasesbyDay2(tt_no_noise[[4]],plottimes,startDate=start_date), f2="0.90")) %>%
  ggplot(aes(x=dates, y=median))+
  geom_line(aes(colour=f2))+
  geom_ribbon(aes(x=dates, ymin=lower25, ymax=upper75, fill=f2), alpha=0.3)+
  scale_fill_manual(values=c("0.65"="#7ec5a3", "0.90"="#f7bb82", "0.50"="#9caee4", "baseline"="#cdcdcd"))+
  scale_colour_manual(values=c("0.65"="#7ec5a3", "0.90"="#f7bb82", "0.50"="#9caee4", "baseline"="#cdcdcd"))+
  geom_vline(xintercept=start_date+firstdiv$results$time, colour=c("#9caee4", "#7ec5a3", "#f7bb82"), size=1, linetype="dashed")+
  geom_vline(xintercept=dmy("17-05-2020"), colour="#cdcdcd", size=1)+
  xlim(dmy("01-05-2020"), dmy("01-07-2020"))+
  ylim(0,600) +
  labs(x="Date", y="Active cases", fill=expression(f[2]), colour=expression(f[2]))+
  theme_minimal()

 p = last_plot()
p
ggsave('Figure2b.pdf', width=6, height=2, units = 'in', dpi = 900)
