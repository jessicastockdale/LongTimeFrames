#----------------------------------------------------------
# preliminaries

rm(list=ls())
library(tidyverse)
library(lubridate)
source('../sd_funcs.R')
library(deSolve)
library(parallel)
set.seed(684645)

if(file.exists('Figure4.RData')){load('Figure4.RData')}else{
#----------------------------------------------------------
# Read in the files
bcdata = read.csv("../bc_casecounts2204.csv", header=TRUE)
bcdata[,1] = dmy(bcdata[,1])
bcdata$Diffs = c(bcdata$Cases[2] - bcdata$Cases[1],
                 diff(bcdata$Cases))
bcdata$day = 1:nrow(bcdata) 


#----------------------------------------------------------
# Various parameters

N = 5.1e6  #population of BC

# initial state of the system
pars = list(N = N, D = 5, R0 = 2.5, k1 = 0.2 ,k2 = 1, q = 0.05, r = 1, ur = 0.4,f = 0.75)
fsi = with(pars,
           r/(r+ur))
nsi = 1 - fsi
i0 = 40 #  we need an initial number of cases, split among the different categories, to set things off

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

# time period
times = seq(from=-20,
            to=25,
            by=0.1)

# timing of social distancing
mysdtiming=function(t, pars){
  ifelse( (t > startTime & t< stopTime),
          pars$f,
          1)
}

# R0 = 2.57 from pre-distancing model fit.
# and psi = 1.94

#-----------------------------------------------------------

Nreps = 50

R0_vec = rnorm(Nreps, mean=2.57, sd=0.05) # sd 0.05 for distancing-on version

times = c(seq(-20.0, max(bcdata$day), 0.1))

startTime = 18 # start of distancing
stopTime = 1000

pars = list(N = N, D = 5, R0 = 2.57, k1 = 0.2 ,k2 = 1, q = 0.05, r = 1, ur = 0.4,f = 0.75)

fit_f = function(state,
                 pars,
                 R0,
                 data,
                 ratio){
  library(deSolve)
  pars$R0 = R0
  MLEs = rep(1, dim(data)[1])
  for( j in 14:dim(data)[1] ){
    ffit = optimize(function(f){ negloglikef(f,
                                             ratio,
                                             pars,
                                             data[1:j,],
                                             state,
                                             times[1:which(times==j)], mysdtiming, overdisp = 5.0) }, lower=0, upper=1)
    MLEs[j] = ffit$minimum
  }
  
  return( MLEs )
}

ncore = 4 # user would need to change based on their #cores
cl = makeCluster( ncore )
clusterExport(cl, varlist=ls(), envir=environment())

system.time({
results = parLapply( cl, X=1:Nreps, function(i){fit_f(state,
                                       pars,
                                       R0_vec[i],
                                       bcdata,
                                       1.94)})}) # ratio = 1.94
stopCluster(cl)

t = dim(bcdata)[1]
res = cbind(1:t, results[[1]])
for( i in 2:Nreps )
  res = rbind(res, cbind(1:t,results[[i]]))

res = data.frame(res)
names(res) = c('t','mle')
temp = res %>% group_by(t) %>% summarise(mle=median(mle))

# We want to 'accept' the MLE when it changes by less than 5% over 3 successive days
for (i in 14:dim(bcdata)[1]){
  if (temp$mle[i]<0.999 && temp$mle[i]>0.0001){ # to stop it accepting right away - need a better fix than this
    if (abs((temp$mle[i+1]-temp$mle[i])/temp$mle[i])<0.05){ # compare ith to the i+1th
      if (abs((temp$mle[i+2]-temp$mle[i+1])/temp$mle[i+1])<0.05){
        # accept that set of 3 days
        return()
      }
    }
  }
}
# the result:
i
rec_start = i
rec_end = i+2

save.image("Figure4.RData")
}

 plot_data = res %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975))
 
 
 p1 <- ggplot(data = plot_data, aes(x=t,y=median))+
  geom_ribbon(data = plot_data[1:rec_start,], aes(x=t[1:rec_start], ymin=L[1:rec_start], ymax=U[1:rec_start]), alpha=0.5, fill = "gray50") +
  geom_ribbon(data = plot_data[-(1:rec_start-1),], aes(x=t, ymin=L, ymax=U), alpha=0.5, fill = "skyblue1") +
  
  geom_point(data = plot_data[1:rec_start,], aes(x=t, y=median),size=1.6,shape=17, color="gray39") +
  geom_point(data = plot_data[-(1:rec_start-1),], aes(x=t, y=median),size=1.6,shape=17, color="royalblue3") +
    
  labs(x='Time (days)', y="f (MLE)") +
  theme_minimal() +
    
  geom_vline(xintercept = as.numeric(res$t[18]), color = "gray15", size=1.5, alpha = 0.5) +
  annotate(geom = "text",
           x = res$t[17]-0.75,
           y = min(res$mle)+0.6,
           label = "Physical distancing",
           color = "black", angle = 90) +
  annotate("rect", xmin=as.numeric(res$t[rec_start]), xmax = as.numeric(res$t[rec_end]), 
           ymin=0,ymax=1 ,fill="royalblue3", alpha = 0.5) +
  annotate(geom = "text", x = res$t[rec_start-1] - 0.75, y = 0.75,
           label = "MLE accepted", color="royalblue3", angle = 90)  +
  theme_minimal()
res
#ggsave('Figure4.pdf', width=6, height=3, units = 'in', dpi = 900)

p2 <- bcdata %>%
  mutate(Time = as.numeric(Date)-min(as.numeric(Date))) %>%
  ggplot(aes(x=Time,y=Diffs))+geom_point()+geom_line() + 
  labs(x='', y="Daily cases") +theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#ggsave('Figure4_casecounts.pdf', width=6, height=3, units = 'in', dpi = 900)


library(gridExtra )
grid.arrange(grobs = list(p2, p1), nrow = 2, heights = c(1,3))
# saved as 6x4 inches
