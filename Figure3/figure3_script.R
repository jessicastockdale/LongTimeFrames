rm(list=ls())
library(deSolve)
library(tidyverse)
library(lubridate)
library(reshape2)
library(parallel)
library(nonpar)
library(ggpubr)

source("sd_funcsfig3.R") # has appropriate startdate, quantiles and sd
set.seed(1211)


if(file.exists("figure3_output.RData")==FALSE){
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
i0 = 470 #  we need an initial number of cases, split among the different categories E1, E2 I, to set things off - we take the number of cases identified in BC on the 11 days preceeding April 10th (average length of time in E1+E2+I)

state = c(S = nsi*(N-i0), # We don't reduce S, since an insignificant no. of ppl have been infected
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



# social distancing first implemented on April 4th
start_date = dmy("10-04-2020")

# social distancing relaxed starting May 17th
relax_date = dmy("17-05-2020")

# number of days between social distancing implemented and relaxed
t_relax = as.numeric(relax_date)-as.numeric(start_date)

## for line plot ##

ss = c(pars$f, 0.6, 0.7, 0.8, 0.9, 1.0)

niter = 100

nReps = 100

ncore = detectCores()
cl = makeCluster( ncore )

fd1=fd2=fd3=fd4=list()
dd = list()
dd[[1]] = dd[[2]] = dd[[3]] = dd[[4]] = matrix(NA, nrow=length(ss)-1, ncol = niter)

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
                                                                 timefunc=timelist2[[x]]))
  fd1[[j]] = firstdiverge2(tt, increasing=T, cc=5)
  fd2[[j]] = firstdiverge2(tt, increasing=T, cc=10)
  fd3[[j]] = firstdiverge2(tt, increasing=T, cc=15)
  fd4[[j]] = firstdiverge2(tt, increasing=T, cc=20)

  dd[[1]][,j] = fd1[[j]]$results$time - t_relax
  dd[[2]][,j] = fd2[[j]]$results$time - t_relax
  dd[[3]][,j] = fd3[[j]]$results$time - t_relax
  dd[[4]][,j] = fd4[[j]]$results$time - t_relax
  }
stopCluster(cl)

save.image(file="figure3_output.RData")
} else(load("figure3_output.RData"))

load("figure3_output.RData")

  qq = list()
  qq[[1]] = qq[[2]] = qq[[3]] = qq[[4]] = matrix(NA, nrow= nrow(dd[[1]]), ncol=3)


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

  # Violin plot

  fs = c(0.6,0.7,0.8,0.9,1)
  days5 = data.frame(days=matrix(t(dd[[1]]), ncol=1, byrow=T), f=rep(fs, each=100), cc=rep("5",500)) # days with cc=5, f2 = c(1,0.9,0.8,0.7,0.6)
  days10 = data.frame(days=matrix(t(dd[[2]]), ncol=1, byrow=T), f=rep(fs, each=100), cc=rep("10",500)) # days with cc=10, f2 = c(1,0.9,0.8,0.7,0.6)
  days15 = data.frame(days=matrix(t(dd[[3]]), ncol=1, byrow=T), f=rep(fs, each=100), cc=rep("15",500)) # days with cc=15, f2 = c(1,0.9,0.8,0.7,0.6)
  days20 = data.frame(days=matrix(t(dd[[4]]), ncol=1, byrow=T), f=rep(fs, each=100), cc=rep("20",500)) # days with cc=20, f2 = c(1,0.9,0.8,0.7,0.6)

  df = rbind(days5, days10, days15, days20)
  df$days[is.na(df$days)==TRUE] = max(df$days[is.na(df$days)==FALSE])

ccs = unique(df$cc)
Threshold = factor(df$cc, levels = ccs)
  g6 = ggplot(df)+
    geom_violin(aes(x = factor(f, levels=c("1","0.9","0.8","0.7","0.6")), y=days,
                    fill=Threshold, colour=Threshold))+
    labs(x= c("f2"), y="Days until threshold")+ theme_minimal() +
    ylim(c(0,80)) +
      scale_colour_manual(values=c("#7ec5a3", "#f7bb82", "#9caee4", "#cdcdcd")) +
     scale_fill_manual(values=c("#7ec5a3", "#f7bb82", "#9caee4", "#cdcdcd"))

      g6


 ggsave('Figure3.pdf', width=6, height=3, units = 'in', dpi = 900)
