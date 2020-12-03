rm(list=ls())
library(deSolve)
library(tidyverse)
library(lubridate)
library(reshape2)
library(parallel)
library(nonpar)
source("sd_funcs.R") # Version has correct sd, quantiles etc.
set.seed(1211)

if ( file.exists('FigureS1T.RData') ){load('FigureS1T.RData')}else{
  sds = seq(0.01, 0.4, by = 0.05)
  zz = list()
  yy = list()

  ncore = 32
  cl = makeCluster(ncore)
  for (ii in 1:length(sds)) {
    cat(date())
    cat(sprintf(' %d / %d\n', ii, length(sds)))
    N = 5.1e6  #population of BC
    sd = sds[ii]
    pars = list(N = N,
                D = 5,
                R0 = 2.57,
                k1 = 1/4,
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

    nReps = 1000

    time_list = list()
    longtimes = seq(from = -20,
                    to = 100,
                    by = 0.1) # how long to simulate for


    # general function for social distancing timing
    timing = function( t, f1, f2, tstar ){
      ifelse( t<tstar, f1, f2)
    }


    # social distancing first implemented on April 4th
    start_date = dmy("10-04-2020")

    # social distancing relaxed starting May 17th
    relax_date = dmy("17-05-2020")

    # number of days between social distancing implemented and relaxed
    t_relax = as.numeric(relax_date)-as.numeric(start_date)


    # No relaxation of distancing : f stays constant at 0.36
    time_list[[1]] = function(t, pars) {timing(t, pars$f, pars$f, t_relax)} # always 0.36

    # Weak relaxation : f -> 0.5 at relax_date
    time_list[[2]] = function(t, pars) {timing(t, pars$f, 0.6, t_relax)}

    # Medium relaxation : f -> 0.65 at relax_date
    time_list[[3]] = function(t, pars) {timing(t, pars$f, 0.7, t_relax)}

    # Strong relaxation : f => 0.9 at relax_date
    time_list[[4]] = function(t, pars) {timing(t, pars$f, 0.8, t_relax)}

    # Strong relaxation : f => 0.9 at relax_date
    time_list[[5]] = function(t, pars) {timing(t, pars$f, 0.9, t_relax)}

    clusterExport(cl, varlist=ls())
    tt2 = parLapply(cl, X=time_list, function(x) multisolve(params=pars,
                                                            state,
                                                            longtimes,
                                                            nReps = nReps,
                                                            timefunc=x, sd = sd))




    names(tt2) = c(0.36, 0.6, 0.7, 0.8, 0.9)

    yy[[ii]] = firstdiverge2(tt2, increasing=T, cc =  10)
    rm(list=c('tt2'))
  }

  stopCluster(cl)

  save.image(file='FigureS1T.RData')
}

dd = matrix(NA, nrow=length(sds), ncol=4)
for (i in 1:length(yy)){
  dd[i,] = (yy[[i]]$results$time-t_relax)
}

ys = data.frame(matrix(dd, ncol=1)) %>%
  mutate(f = c(rep("0.6", length(sds)), rep("0.7", length(sds)), rep("0.8", length(sds)),  rep("0.9", length(sds)))) %>%
  mutate(sd = rep(sds, 4))
colnames(ys) = c("days", "f2", "sd")

vplot = ggplot(data = ys)+
  geom_line(aes(x=sd, y=days, colour=f2, linetype = f2), size=0.5)+
geom_point(aes(x=sd, y=days, colour=f2, shape = f2), size=1)

vplot + xlim(c(0, 0.4)) + scale_y_continuous(limits=c(0, 60)) + theme_minimal() + ylab('Days until threshold') + xlab('sd')

vplot = last_plot()
vplot
ggsave('FigureS1top.pdf', width=6, height=3, units = 'in', dpi = 900)
