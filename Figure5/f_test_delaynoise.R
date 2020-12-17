#----------------------------------------------------------
# preliminaries

rm(list=ls())
library(tidyverse)
library(lubridate)
source('../sd_funcs.R')
library(deSolve)
library(parallel)
set.seed(684645)

if(file.exists('noisedelay.RData')){load('noisedelay.RData')}
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
# initial pars
pars = list(N = N,
            D = 5,
            R0 = 2.5,
            k1 = 1/5,
            k2 = 1,
            q = 0.05,
            r = 1,
            ur = 0.4,
            f = 0.75)

# initial state of the system
fsi = with(pars,
           r/(r+ur))
nsi = 1 - fsi
i0 = 40 
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

# Changing level of noise:

Nreps = 50

R0_vec = rnorm(Nreps, mean=2.57, sd=0.05)

# simulate the time series:
opt_pars = list(N = N, D = 5, R0 = 2.57, k1 = 1/5, k2 = 1, q = 0.05, r = 1, ur = 0.4, f = 0.36)
# We use our most recent f MLE from the last part
times = c(seq(-20.0, max(bcdata$day), 0.1))
startTime = 18 
stopTime = 1000

out = as.data.frame(ode(y = state,times = times,func = socdistmodel,parms=opt_pars,sdtiming = mysdtiming))

# We need the same format as 'bcdata' - 4 columns Date, Cases, Diffs, day
simdata = bcdata
# For simplicity, we initialise the simulation to 8 cases as observed in reality by day 0
simdata$Diffs[1] = round(getlambd(out,opt_pars,day=1,ratio = 1.94, data = bcdata,sampFrac = 0.35,
                                  delayShape = 1.73, delayScale = 9.85, change_days = "2020-03-14"))
simdata$Cases[1] = 8 + simdata$Diffs[1]

for (i in 2:length(simdata$day)){
  # gives the model number of observed cases on a given day e.g. our simulated cases per day
  obs = getlambd(out,opt_pars,day=i,ratio = 1.94, data = bcdata,sampFrac = 0.35,
                 delayShape = 1.73, delayScale = 9.85, change_days = "2020-03-14") 
  # we round to the nearest whole number of cases
  simdata$Cases[i] = simdata$Cases[i-1] + round(obs)
  simdata$Diffs[i] = round(obs)
}

# Second noise level to test:
simdata2 = bcdata
phi = 5
dScale = 9.85
dShape = 1.73

simdata2$Diffs[1] = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=1,ratio = 1.94, data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
simdata2$Cases[1] = 8 + simdata2$Diffs[1]
for (i in 2:length(simdata2$day)){
  # gives the model number of observed cases on a given day e.g. our simulated cases per day
  obs = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=i,ratio = 1.94, data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
  # we need to round the nearest whole number of cases
  simdata2$Cases[i] = simdata2$Cases[i-1] + round(obs)
  simdata2$Diffs[i] = round(obs)
}


# Third noise level to test:
simdata3 = bcdata
phi = 10
dScale = 9.85
dShape = 1.73

simdata3$Diffs[1] = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=1,ratio = 1.94, data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
simdata3$Cases[1] = 8 + simdata3$Diffs[1]
for (i in 2:length(simdata3$day)){
  # gives the model number of observed cases on a given day e.g. our simulated cases per day
  obs = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=i,ratio = 1.94, data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
  # we need to round the nearest whole number of cases
  simdata3$Cases[i] = simdata3$Cases[i-1] + round(obs)
  simdata3$Diffs[i] = round(obs)
}



 fit_f = function(state,
                 pars,
                 R0,
                 data1, data2, data3,
                 ratio,
                 overdisp1, overdisp2, overdisp3){
   # data1 and data2 need to be the same size
  library(deSolve)
  pars$R0 = R0
  MLEs1 = rep(1, dim(data1)[1])
  MLEs2 = rep(1, dim(data2)[1])
  MLEs3 = rep(1, dim(data3)[1])
  for( j in 14:dim(data1)[1] ){
        ffit=optimize(function(f) negloglikef(f,
                                              ratio,
                                              pars,
                                              data1[1:j,],
                                              state,
                                              times[1:which(times==j)],
                                              mysdtiming,  overdisp = overdisp1),
                          lower=0.0, upper=1.0)
        ffit2=optimize(function(f) negloglikef(f,
                                              ratio,
                                              pars,
                                              data2[1:j,],
                                              state,
                                              times[1:which(times==j)],
                                              mysdtiming,  overdisp = overdisp2),
                      lower=0.0, upper=1.0)
        ffit3=optimize(function(f) negloglikef(f,
                                               ratio,
                                               pars,
                                               data3[1:j,],
                                               state,
                                               times[1:which(times==j)],
                                               mysdtiming,  overdisp = overdisp3),
                       lower=0.0, upper=1.0)
    MLEs1[j] = ffit$minimum
    MLEs2[j] = ffit2$minimum
    MLEs3[j] = ffit3$minimum
  }
  
  return( cbind(MLEs1, MLEs2, MLEs3) )
}

ncore = 4
cl = makeCluster( ncore )
clusterExport(cl, varlist=ls(), envir=environment())

system.time({
results = parLapply( cl, X=1:Nreps,
                     function(i){
                       fit_f(state = state,
                       pars = opt_pars,
                       R0 = R0_vec[i],
                       data1 = simdata, data2 = simdata2, data3 = simdata3,
                       ratio = 1.94,
                       overdisp1 = 50, overdisp2=5, overdisp3 = 10)})})
stopCluster(cl)

t = dim(simdata)[1]

res = cbind(1:t, results[[1]][,1])
for( i in 2:Nreps )
  res = rbind(res, cbind(1:t,results[[i]][,1]))

res = data.frame(res)
names(res) = c('t','mle')
temp1 = res %>% group_by(t) %>% summarise(mle=median(mle))


res2 = cbind(1:t, results[[1]][,2])
for( i in 2:Nreps )
  res2 = rbind(res2, cbind(1:t,results[[i]][,2]))

res2 = data.frame(res2)
names(res2) = c('t','mle')
temp2 = res2 %>% group_by(t) %>% summarise(mle=median(mle))

res3 = cbind(1:t, results[[1]][,3])
for( i in 2:Nreps )
  res3 = rbind(res3, cbind(1:t,results[[i]][,3]))

res3 = data.frame(res3)
names(res3) = c('t','mle')
temp3 = res3 %>% group_by(t) %>% summarise(mle=median(mle))



# We want to 'accept' the MLE when it changes by less than 5% over 3 successive days
for (i in 14:dim(simdata)[1]){
  if (temp1$mle[i]<0.999 && temp1$mle[i]>0.0001){ # to stop it accepting right away - need a better fix than this
    if (abs((temp1$mle[i+1]-temp1$mle[i])/temp1$mle[i])<0.05){ # compare ith to the i+1th
      if (abs((temp1$mle[i+2]-temp1$mle[i+1])/temp1$mle[i+1])<0.05){
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

for (i in 14:dim(simdata2)[1]){
  if (temp2$mle[i]<0.999 && temp2$mle[i]>0.0001){ # to stop it accepting right away - need a better fix than this
    if (abs((temp2$mle[i+1]-temp2$mle[i])/temp2$mle[i])<0.05){ # compare ith to the i+1th
      if (abs((temp2$mle[i+2]-temp2$mle[i+1])/temp2$mle[i+1])<0.05){
        # accept that set of 3 days
        return()
      }
    }
  }
}
# the result:
i
rec_start2 = i
rec_end2 = i+2

for (i in 14:dim(simdata3)[1]){
  if (temp3$mle[i]<0.999 && temp3$mle[i]>0.0001){ # to stop it accepting right away - need a better fix than this
    if (abs((temp3$mle[i+1]-temp3$mle[i])/temp3$mle[i])<0.05){ # compare ith to the i+1th
      if (abs((temp3$mle[i+2]-temp3$mle[i+1])/temp3$mle[i+1])<0.05){
        # accept that set of 3 days
        return()
      }
    }
  }
}
# the result:
i
rec_start3 = i
rec_end3 = i+2


save.image("noisedelay.RData")

# rename to avoid ggplot error later...
res2a <- res2
res3a <- res3
resa <- res
rec_starta <- rec_start
rec_enda <- rec_end
rec_start3a <- rec_start3
rec_end3a <- rec_end3
rec_start2a <- rec_start2
rec_end2a <- rec_end2

tib2a = res2a %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975))

tib3a = res3a %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975))


tiba = resa %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975)) 


  p3 <- ggplot(data = tiba, aes(x=t,y=median))+
  geom_ribbon(data = tiba[1:rec_starta,], aes(x=t[1:rec_starta], ymin=L[1:rec_starta], ymax=U[1:rec_starta]), alpha=0.5, fill="gray") +
  geom_ribbon(data=tib3a[1:rec_start3a,], aes(x=tib3a$t[1:rec_start3a], ymin=tib3a$L[1:rec_start3a], ymax=tib3a$U[1:rec_start3a]), alpha=0.5, fill="gray") +
  geom_ribbon(data=tib2a[1:rec_start2a,], aes(x=tib2a$t[1:rec_start2a], ymin=tib2a$L[1:rec_start2a], ymax=tib2a$U[1:rec_start2a]), alpha=0.5, fill="gray") +
    
  geom_ribbon(data = tiba[-(1:rec_starta-1),], aes(x=tiba$t[-(1:rec_starta-1)], ymin=tiba$L[-(1:rec_starta-1)], ymax=tiba$U[-(1:rec_starta-1)]), alpha=0.5, fill="palegreen1") +
  geom_ribbon(data = tib3a[-(1:rec_start3a-1),], aes(x=tib3a$t[-(1:rec_start3a-1)], ymin=tib3a$L[-(1:rec_start3a-1)], ymax=tib3a$U[-(1:rec_start3a-1)]), alpha=0.5, fill="gold1") +
  geom_ribbon(data=tib2a[-(1:rec_start2a-1),], aes(x=tib2a$t[-(1:rec_start2a-1)], ymin=tib2a$L[-(1:rec_start2a-1)], ymax=tib2a$U[-(1:rec_start2a-1)]), alpha=0.5, fill="skyblue1") +
    
  geom_point(data = tiba[1:rec_starta,], size=1.8,shape=18, color="gray39") +
  geom_point(data=tib3a[1:rec_start3a,], x=tib3a$t[1:rec_start3a],y=tib3a$median[1:rec_start3a],size=1.6,shape=19, color="gray39") +
  geom_point(data=tib2a[1:rec_start2a,], x=tib2a$t[1:rec_start2a],y=tib2a$median[1:rec_start2a],size=1.6,shape=17, color="gray39") +
  geom_point(data=tib2a[-(1:rec_start2a-1),], x=tib2a$t[-(1:rec_start2a-1)],y=tib2a$median[-(1:rec_start2a-1)],size=1.6,shape=17, color="royalblue3") +
  geom_point(data=tib3a[-(1:rec_start3a-1),], x=tib3a$t[-(1:rec_start3a-1)],y=tib3a$median[-(1:rec_start3a-1)],size=1.6,shape=19, color="darkorange2") +
  geom_point(data = tiba[-(1:rec_starta-1),], size=1.8,shape=18, color="springgreen4") +
     
  labs(x='Time (days)', y="f (MLE)") +
  theme_minimal() +
  geom_vline(xintercept = as.numeric(resa$t[18]), color = "gray25", size=1.5, alpha = 0.5) +
  annotate(geom = "text",
           x = resa$t[17]-0.75,
           y = min(resa$mle)+0.6,
           label = "Physical distancing",
           color = "black", angle = 90) +
  
   annotate("rect", xmin=as.numeric(resa$t[rec_starta]), xmax = as.numeric(resa$t[rec_enda]), 
            ymin=0,ymax=1 ,fill="springgreen4", alpha = 0.5) +
   annotate("rect", xmin=as.numeric(res2a$t[rec_start2a]), xmax = as.numeric(res2a$t[rec_end2a]), 
             ymin=0,ymax=1 ,fill="royalblue3", alpha = 0.5) +
   annotate("rect", xmin=as.numeric(res3a$t[rec_start3a]), xmax = as.numeric(res3a$t[rec_end3a]), 
             ymin=0,ymax=1 ,fill="darkorange2", alpha = 0.5) +
    
  annotate(geom = "text", x = resa$t[rec_starta-1] - 0.75, y = 0.75,
           label = "MLE accepted", color = "springgreen4", angle = 90)  +
  annotate(geom = "text", x = res2a$t[rec_start2a-1] - 0.75, y = 0.75,
           label = "MLE accepted", color = "royalblue3", angle = 90) +
  annotate(geom = "text", x = res3a$t[rec_start3a-1] - 0.75, y = 0.75,
             label = "MLE accepted", color = "darkorange2", angle = 90)+ 
    theme(plot.margin = unit(c(0.2,0,1,0), "cm"))
  # (add margin for grid arrange)
  
  case_data <- simdata %>%
    mutate(Time = as.numeric(Date)-min(as.numeric(Date)))
  case_data <- case_data[,c(5,3)]
  case_data <- cbind(case_data, simdata2[,3], simdata3[,3])
  names(case_data) = c("Time", "Diffs1", "Diffs2", "Diffs3")
  
  p4 <- ggplot(data = case_data, aes(x=Time,y=Diffs1))+
    geom_point(aes(x=Time,y=Diffs2), color = "royalblue3")+geom_line(aes(x=Time,y=Diffs2), color = "royalblue3") + 
    geom_point(aes(x=Time,y=Diffs3),color = "darkorange2")+geom_line(aes(x=Time,y=Diffs3),color = "darkorange2") + 
    geom_point(color = "springgreen4")+geom_line(color = "springgreen4") + 
    labs(x='', y="Daily cases") +theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))
  # (add margin for grid arrange)
  
  
  grid.arrange(grobs = list(p4, p3), nrow = 2, heights = c(1,3))
  # saved as 6x4inches 
  


if(file.exists('delay.RData')){load('delay.RData')}
#-----------------------------------------------------------

# Changing length of delay:
set.seed(171094)
  
Nreps = 50

R0_vec = rnorm(Nreps, mean=2.57, sd=0.05)

# simulate the time series:
opt_pars = list(N = N, D = 5, R0 = 2.57, k1 = 1/5, k2 = 1, q = 0.05, r = 1, ur = 0.4, f = 0.36)
# We use our most recent f MLE from the last part
times = c(seq(-20.0, max(bcdata$day), 0.1))
startTime = 18 
stopTime = 1000

# First delay length to test:
simdata = bcdata
# Let's keep the noise fixed to phi = 5
phi = 5
dScale = 9.85
dShape = 1.73
meanDelay = dScale * gamma(1 + 1 / dShape) # = 
meanDelay
varDelay = dScale^2*(gamma(1 + 2/dShape) - (gamma(1 + 1/dShape))^2)
varDelay

out = as.data.frame(ode(y = state,times = times,func = socdistmodel,parms=opt_pars,sdtiming = mysdtiming))

simdata$Diffs[1] = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=1,ratio = 1.94,data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
simdata$Cases[1] = 8 + simdata$Diffs[1]
for (i in 2:length(simdata$day)){
  # gives the model number of observed cases on a given day e.g. our simulated cases per day
  obs = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=i,ratio = 1.94,data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
  # we need to round the nearest whole number of cases
  simdata$Cases[i] = simdata$Cases[i-1] + round(obs)
  simdata$Diffs[i] = round(obs)
}

# Second delay length to test, around 4 days instead of 8:
simdata2 = bcdata
phi = 5
dScale = 5
dShape = 2 # mean ~4.4, var ~5
meanDelay = dScale * gamma(1 + 1 / dShape) # = 
varDelay = dScale^2*(gamma(1 + 2/dShape) - (gamma(1 + 1/dShape))^2)

simdata2$Diffs[1] = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=1,ratio = 1.94,data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
simdata2$Cases[1] = 8 + simdata2$Diffs[1]
for (i in 2:length(simdata2$day)){
  # gives the model number of observed cases on a given day e.g. our simulated cases per day
  obs = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=i,ratio = 1.94,data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
  # we need to round the nearest whole number of cases
  simdata2$Cases[i] = simdata2$Cases[i-1] + round(obs)
  simdata2$Diffs[i] = round(obs)
}

# Third delay length to test (very small delay):
simdata3 = bcdata
phi = 5
dScale = 1
dShape = 10
meanDelay = dScale * gamma(1 + 1 / dShape) # = 
meanDelay
varDelay = dScale^2*(gamma(1 + 2/dShape) - (gamma(1 + 1/dShape))^2)
varDelay

simdata3$Diffs[1] = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=1,ratio = 1.94,data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
simdata3$Cases[1] = 8 + simdata3$Diffs[1]
for (i in 2:length(simdata3$day)){
  # gives the model number of observed cases on a given day e.g. our simulated cases per day
  obs = round(rnbinom(1, size=phi, mu=getlambd(out,opt_pars,day=i,ratio = 1.94,data = bcdata,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
  # we need to round the nearest whole number of cases
  simdata3$Cases[i] = simdata3$Cases[i-1] + round(obs)
  simdata3$Diffs[i] = round(obs)
}



fit_f = function(state,
                 pars,
                 R0,
                 data1, data2, data3,
                 ratio,
                 overdisp1, overdisp2, overdisp3, dShape1, dScale1, dShape2, dScale2, dShape3, dScale3){
  # data1 and data2 need to be the same size
  library(deSolve)
  pars$R0 = R0
  MLEs1 = rep(1, dim(data1)[1])
  MLEs2 = rep(1, dim(data2)[1])
  MLEs3 = rep(1, dim(data3)[1])
  for( j in 14:dim(data1)[1] ){
    ffit=optimize(function(f) negloglikef(f,
                                          ratio,
                                          pars,
                                          data1[1:j,],
                                          state,
                                          times[1:which(times==j)],
                                          mysdtiming,  overdisp = overdisp1,
                                          delayShape = dShape1, delayScale = dScale1),
                  lower=0.0, upper=1.0)
    ffit2=optimize(function(f) negloglikef(f,
                                           ratio,
                                           pars,
                                           data2[1:j,],
                                           state,
                                           times[1:which(times==j)],
                                           mysdtiming,  overdisp = overdisp2,
                                           delayShape = dShape2, delayScale = dScale2),
                   lower=0.0, upper=1.0)
    ffit3=optimize(function(f) negloglikef(f,
                                           ratio,
                                           pars,
                                           data3[1:j,],
                                           state,
                                           times[1:which(times==j)],
                                           mysdtiming,  overdisp = overdisp3,
                                           delayShape = dShape3, delayScale = dScale3),
                   lower=0.0, upper=1.0)
    MLEs1[j] = ffit$minimum
    MLEs2[j] = ffit2$minimum
    MLEs3[j] = ffit3$minimum
  }
  
  return( cbind(MLEs1, MLEs2, MLEs3) )
}

ncore = 4
cl = makeCluster( ncore )
clusterExport(cl, varlist=ls(), envir=environment())

phi
system.time({
  results = parLapply( cl, X=1:Nreps, function(i){fit_f(state = state,
                                                        pars = opt_pars,
                                                        R0 = R0_vec[i],
                                                        data1 = simdata, data2 = simdata2, data3 = simdata3,
                                                        ratio = myfit$par[2],
                                                        overdisp1 = phi, overdisp2 = phi, overdisp3 = phi,
                                                        dShape1 = 1.73, dScale1 = 9.85, dShape2 = 2,
                                                        dScale2 = 5, dShape3 = 10, dScale3 = 1)})})
stopCluster(cl)

t = dim(simdata)[1]

res = cbind(1:t, results[[1]][,1])
for( i in 2:Nreps )
  res = rbind(res, cbind(1:t,results[[i]][,1]))

res = data.frame(res)
names(res) = c('t','mle')
temp1 = res %>% group_by(t) %>% summarise(mle=median(mle))

res2 = cbind(1:t, results[[1]][,2])
for( i in 2:Nreps )
  res2 = rbind(res2, cbind(1:t,results[[i]][,2]))

res2 = data.frame(res2)
names(res2) = c('t','mle')
temp2 = res2 %>% group_by(t) %>% summarise(mle=median(mle))

res3 = cbind(1:t, results[[1]][,3])
for( i in 2:Nreps )
  res3 = rbind(res3, cbind(1:t,results[[i]][,3]))

res3 = data.frame(res3)
names(res3) = c('t','mle')
temp3 = res3 %>% group_by(t) %>% summarise(mle=median(mle))


# We want to 'accept' the MLE when it changes by less than 5% over 3 successive days
for (i in 14:dim(simdata)[1]){
  if (temp1$mle[i]<0.999 && temp1$mle[i]>0.0001){ # to stop it accepting right away 
    if (abs((temp1$mle[i+1]-temp1$mle[i])/temp1$mle[i])<0.05){ # compare ith to the i+1th
      if (abs((temp1$mle[i+2]-temp1$mle[i+1])/temp1$mle[i+1])<0.05){
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

for (i in 14:dim(simdata2)[1]){
  if (temp2$mle[i]<0.999 && temp2$mle[i]>0.0001){ # to stop it accepting right away 
    if (abs((temp2$mle[i+1]-temp2$mle[i])/temp2$mle[i])<0.05){ # compare ith to the i+1th
      if (abs((temp2$mle[i+2]-temp2$mle[i+1])/temp2$mle[i+1])<0.05){
        # accept that set of 3 days
        return()
      }
    }
  }
}
# the result:
i
rec_start2 = i
rec_end2 = i+2


for (i in 14:dim(simdata3)[1]){
  if (temp3$mle[i]<0.999 && temp3$mle[i]>0.0001){ # to stop it accepting right away 
    if (abs((temp3$mle[i+1]-temp3$mle[i])/temp3$mle[i])<0.05){ # compare ith to the i+1th
      if (abs((temp3$mle[i+2]-temp3$mle[i+1])/temp3$mle[i+1])<0.05){
        # accept that set of 3 days
        return()
      }
    }
  }
}
# the result:
i
rec_start3 = i
rec_end3 = i+2


save.image("delay.RData")
 


tib2 = res2 %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975))

tib3 = res3 %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975))

tib = res %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975)) 


p1 <- ggplot(data = tib, aes(x=t,y=median))+
  geom_ribbon(data = tib[1:rec_start,], aes(x=t[1:rec_start], ymin=L[1:rec_start], ymax=U[1:rec_start]), alpha=0.5, fill="gray") +
  geom_ribbon(data=tib3[1:rec_start3,], aes(x=tib3$t[1:rec_start3], ymin=tib3$L[1:rec_start3], ymax=tib3$U[1:rec_start3]), alpha=0.5, fill="gray") +
  geom_ribbon(data=tib2[1:rec_start2,], aes(x=tib2$t[1:rec_start2], ymin=tib2$L[1:rec_start2], ymax=tib2$U[1:rec_start2]), alpha=0.5, fill="gray") +
  
  geom_ribbon(data = tib[-(1:rec_start-1),], aes(x=tib$t[-(1:rec_start-1)], ymin=tib$L[-(1:rec_start-1)], ymax=tib$U[-(1:rec_start-1)]), alpha=0.5, fill="skyblue1") +
  geom_ribbon(data = tib3[-(1:rec_start3-1),], aes(x=tib3$t[-(1:rec_start3-1)], ymin=tib3$L[-(1:rec_start3-1)], ymax=tib3$U[-(1:rec_start3-1)]), alpha=0.5, fill="palegreen1") +
  geom_ribbon(data=tib2[-(1:rec_start2-1),], aes(x=tib2$t[-(1:rec_start2-1)], ymin=tib2$L[-(1:rec_start2-1)], ymax=tib2$U[-(1:rec_start2-1)]), alpha=0.5, fill="gold1") +
  
  geom_point(data = tib[1:rec_start,], size=1.8,shape=18, color="gray39") +
  geom_point(data=tib3[1:rec_start3,], x=tib3$t[1:rec_start3],y=tib3$median[1:rec_start3],size=1.6,shape=19, color="gray39") +
  geom_point(data=tib2[1:rec_start2,], x=tib2$t[1:rec_start2],y=tib2$median[1:rec_start2],size=1.6,shape=17, color="gray39") +
  geom_point(data=tib2[-(1:rec_start2-1),], x=tib2$t[-(1:rec_start2-1)],y=tib2$median[-(1:rec_start2-1)],size=1.6,shape=17, color="darkorange2") +
  geom_point(data=tib3[-(1:rec_start3-1),], x=tib3$t[-(1:rec_start3-1)],y=tib3$median[-(1:rec_start3-1)],size=1.6,shape=19, color="springgreen4") +
  geom_point(data = tib[-(1:rec_start-1),], size=1.8,shape=18, color="royalblue3") +
  
  labs(x='Time (days)', y="f (MLE)") +
  theme_minimal() +
  geom_vline(xintercept = as.numeric(res$t[18]), color = "gray25", size=1.5, alpha = 0.5) +
  annotate(geom = "text",
           x = res$t[17]-0.75,
           y = min(res$mle)+0.6,
           label = "Physical distancing",
           color = "black", angle = 90) +
   
  annotate("rect", xmin=as.numeric(res$t[rec_start]), xmax = as.numeric(res$t[rec_end]), 
           ymin=0,ymax=1 ,fill="royalblue3", alpha = 0.5) +
  annotate("rect", xmin=as.numeric(res2$t[rec_start2]), xmax = as.numeric(res2$t[rec_end2]), 
           ymin=0,ymax=1 ,fill="darkorange2", alpha = 0.5) +
  annotate("rect", xmin=as.numeric(res3$t[rec_start3]), xmax = as.numeric(res3$t[rec_end3]), 
           ymin=0,ymax=1 ,fill="springgreen4", alpha = 0.5) +
  
  annotate(geom = "text", x = res$t[rec_start-1] - 0.75, y = 0.75,
           label = "MLE accepted", color = "royalblue3", angle = 90)  +
  annotate(geom = "text", x = res2$t[rec_start2-1] - 0.3, y = 0.75,
           label = "MLE accepted", color = "darkorange2", angle = 90) +
  annotate(geom = "text", x = res3$t[rec_start3-1] - 0.75, y = 0.75,
           label = "MLE accepted", color = "springgreen4", angle = 90) + 
  theme(plot.margin = unit(c(0.2,0,1,0), "cm"))
# (add margin for grid arrange)


case_data <- simdata %>%
  mutate(Time = as.numeric(Date)-min(as.numeric(Date)))
case_data <- case_data[,c(5,3)]
case_data <- cbind(case_data, simdata2[,3], simdata3[,3])
names(case_data) = c("Time", "Diffs1", "Diffs2", "Diffs3")

p2 <- ggplot(data = case_data, aes(x=Time,y=Diffs1))+
  geom_point(aes(x=Time,y=Diffs2), color = "royalblue3")+geom_line(aes(x=Time,y=Diffs2), color = "royalblue3") + 
  geom_point(aes(x=Time,y=Diffs3),color = "darkorange2")+geom_line(aes(x=Time,y=Diffs3),color = "darkorange2") + 
  geom_point(color = "springgreen4")+geom_line(color = "springgreen4") + 
  labs(x='', y="Daily cases") +theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
# (add margin for grid arrange)


grid.arrange(grobs = list(p2, p1), nrow = 2, heights = c(1,3))
# saved as 6x4inches 





# Combine to one figure, margin space already added
grid.arrange(grobs = list(p4, p3, p2, p1), nrow = 4, heights = c(1,3, 1, 3))

