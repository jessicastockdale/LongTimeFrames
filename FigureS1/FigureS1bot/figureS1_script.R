rm(list=ls())
library(deSolve)
library(tidyverse)
library(lubridate)
library(reshape2)
library(parallel)
library(nonpar)
library(ggpubr)
source("sd_funcs_S1.R") # Set to correct sd, quantiles etc.
set.seed(12234911)

if(file.exists("figureS1_output.RData")==FALSE){
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

# uncertainty in:
# D,  k1,  k2, q, r,  ur

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

longtimes = seq(from = -20,
                to = 120, # changed from 80
                by = 0.1) # how long to simulate for

time_list = list()

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

## for line plot ##

ss = c(pars$f, 0.6, 0.7, 0.8, 0.9, 1.0)

niter = 100

nReps = 100

ncore = 32
cl = makeCluster( ncore )

fd1=fd2=fd3=list()

dd = list()
dd[[1]] = matrix(NA, nrow=length(ss)-1, ncol = niter)
  dd[[2]] = matrix(NA, nrow=length(ss)-1, ncol = niter)
  dd[[3]] = matrix(NA, nrow=length(ss)-1, ncol = niter)

for (j in 1:niter){
  cat(j)
  cat("\n")
  cat(niter)
  cat("\n")
  cat(date())
  cat("\n")
  bpars2 = list()
  timelist2 = list()
  results =list()

  for (i in 1:length(ss)){
    bpars2[[i]] = pars
    bpars2[[i]]$f = ss[i]
    timelist2[[i]] = function(t, bpars) {timing(t, 0.36, bpars$f, t_relax)}
  }
  clusterExport(cl, varlist=ls())
  tt = parLapply(cl, X=1:length(bpars2), function(x) multisolve(params=bpars2[[x]],
                                                                 state,
                                                                 longtimes,nReps = nReps,
                                                                 timefunc=timelist2[[x]], uncertainty = 0.1))
  fd1[[j]] = firstdiverge2(tt, increasing=T)

   tt = parLapply(cl, X=1:length(bpars2), function(x) multisolve(params=bpars2[[x]],
                                                                 state,
                                                                 longtimes,nReps = nReps,
                                                                 timefunc=timelist2[[x]], uncertainty = 0.2))


  fd2[[j]] = firstdiverge2(tt, increasing=T)


   tt = parLapply(cl, X=1:length(bpars2), function(x) multisolve(params=bpars2[[x]],
                                                                 state,
                                                                 longtimes,nReps = nReps,
                                                                 timefunc=timelist2[[x]], uncertainty = 0.3))



  fd3[[j]] = firstdiverge2(tt, increasing=T)


  dd[[1]][,j] = fd1[[j]]$results$time - t_relax
  dd[[2]][,j] = fd2[[j]]$results$time - t_relax
  dd[[3]][,j] = fd3[[j]]$results$time - t_relax

  }
stopCluster(cl)

save.image(file="figureS1_output.RData")
} else(load("figureS1_output.RData"))

load("figureS1_output.RData")

  qq = list()
  qq[[1]] = matrix(NA, nrow= nrow(dd[[1]]), ncol=3)
  qq[[2]] = matrix(NA, nrow= nrow(dd[[1]]), ncol=3)
  qq[[3]] = matrix(NA, nrow= nrow(dd[[1]]), ncol=3)

  dd2 = dd

  for (j in 1:length(dd2)){
    for (i in 1:nrow(dd2[[j]])){
      if(all(is.na(dd2[[j]][i,]))==FALSE & sum(is.na(dd2[[j]][i,])) > ncol(dd2[[j]]) - 3){
        dd2[[j]][i,][is.na(dd2[[j]][i,])==FALSE] = NA
      }
    }
  }


  for (i in 1:length(qq)){
    for (k in 1:nrow(dd2[[1]])){
      qq[[i]][k,] = quantile(dd2[[i]][k,][is.na(dd2[[i]][k,])==FALSE], prob=c(0.05, 0.5, 0.95))
    }
  }


  fs = c(0.6,0.7,0.8,0.9,1)
  days01 = data.frame(days=matrix(t(dd[[1]]), ncol=1, byrow=T), f=rep(fs, each=100), u=rep("0.1",500)) # days with uncertainty=0.1, f2 = c(1,0.9,0.8,0.7,0.6)
  days02 = data.frame(days=matrix(t(dd[[2]]), ncol=1, byrow=T), f=rep(fs, each=100), u=rep("0.2",500)) # days with uncertainty=0.2, f2 = c(1,0.9,0.8,0.7,0.6)
  days03 = data.frame(days=matrix(t(dd[[3]]), ncol=1, byrow=T), f=rep(fs, each=100), u=rep("0.3",500)) # days with uncertainty=0.3, f2 = c(1,0.9,0.8,0.7,0.6)

  df = rbind(days01, days02, days03)
  df$days[is.na(df$days)==TRUE] = max(df$days[is.na(df$days)==FALSE])

  uus = unique(df$u)
  Uncertainty = factor(df$u, levels = uus)
  g6 = ggplot(df)+
    geom_violin(aes(x = factor(f, levels=c("1","0.9","0.8","0.7","0.6")), y=days,
                    color=Uncertainty, fill=Uncertainty))+
    labs(x= c("f2"), y="Days until threshold")+ theme_minimal() +
    ylim(c(0,80)) +
    scale_colour_manual(values=c("#7ec5a3", "#f7bb82", "#9caee4", "#cdcdcd")) +
    scale_fill_manual(values=c("#7ec5a3", "#f7bb82", "#9caee4", "#cdcdcd"))
  g6
  ggsave('FigureS1bot.pdf', width=6, height=3, units = 'in', dpi = 900)
