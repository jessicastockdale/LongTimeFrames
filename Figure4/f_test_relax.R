#----------------------------------------------------------
# preliminaries

rm(list=ls())
library(tidyverse)
library(lubridate)
source('../sd_funcs.R')
library(deSolve)
library(parallel)
library(gridExtra)
set.seed(684645)

if(file.exists('Figure4_relax.RData')){load('Figure4_relax.RData')}
#----------------------------------------------------------
# Read in the files
  bcdata = read.csv("../bc_casecounts2204.csv", header=TRUE)
  bcdata[,1] = dmy(bcdata[,1])
  bcdata$Diffs = c(bcdata$Cases[2] - bcdata$Cases[1],
                   diff(bcdata$Cases))
  bcdata$day = 1:nrow(bcdata) 
  
  
  #----------------------------------------------------------
  # Various parameters
  
  # DE parameters
  N = 5.1e6  #population of BC
  
  # initial state of the system
  pars = list(N = N, D = 5, R0 = 2.5, k1 = 1/5, k2 = 1, q = 0.05, r = 1, ur = 0.4, f = 0.36)
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
  
  # R0 = 2.57 from pre-distancing model fit.
  # and psi = 1.94
  
  #-----------------------------------------------------------
  # Simulate relaxation of distancing
  
  # simulated up to May 1st
  times = seq(from=-20,
              to=61,
              by=0.1)
  opt_pars = list(N = N, D = 5, R0 = 2.57, k1 = 1/5, k2 = 1, q = 0.05, r = 1, ur = 0.4, f = 0.36)
  
  mysdtiming=function(t, pars){
    ifelse( (t > startTime & t< stopTime),
            pars$f,
            0.36)
  }
  startTime = 17 # start of RELAXING distancing - rest of the time it's 0.36
  stopTime = 500
  
  # We will use our most recent f MLE from the last part
  out_b4 = as.data.frame(ode(y = state,times = times,func = socdistmodel,parms=opt_pars,sdtiming = mysdtiming))
  
  # simulate the time series:
  times = c(seq(-20.0, 62, 0.1))
  # We initialise with the state of the model on May 1st
  state = unlist(out_b4[811,-1])
  

  opt_pars1 = list(N = N, D = 5, R0 = 2.57, k1 = 1/5, k2 = 1, q = 0.05, r = 1, ur = 0.4, f = 0.9)
  opt_pars2 = list(N = N, D = 5, R0 = 2.57, k1 = 1/5, k2 = 1, q = 0.05, r = 1, ur = 0.4, f = 0.65)
  opt_pars3 = list(N = N, D = 5, R0 = 2.57, k1 = 1/5, k2 = 1, q = 0.05, r = 1, ur = 0.4, f = 0.5)
  
  out1 = as.data.frame(ode(y = state,times = times,func = socdistmodel,parms=opt_pars1,sdtiming = mysdtiming))
  out2 = as.data.frame(ode(y = state,times = times,func = socdistmodel,parms=opt_pars2,sdtiming = mysdtiming))
  out3 = as.data.frame(ode(y = state,times = times,func = socdistmodel,parms=opt_pars3,sdtiming = mysdtiming))

  
  # We need the same format as 'bcdata' - 4 columns Date, Cases, Diffs, day
  simdata1 = data.frame(matrix(NA, nrow = 62, ncol = 4))
  names(simdata1) = c("Date", "Cases", "Diffs", "day")
  simdata1$day = 1:62
  simdata1$Date = seq(as.Date("2020-05-01"), by = "day", length.out = 62)
  simdata2 = data.frame(matrix(NA, nrow = 62, ncol = 4))
  names(simdata2) = c("Date", "Cases", "Diffs", "day")
  simdata2$day = 1:62
  simdata2$Date = seq(as.Date("2020-05-01"), by = "day", length.out = 62)
  simdata3 = data.frame(matrix(NA, nrow = 62, ncol = 4))
  names(simdata3) = c("Date", "Cases", "Diffs", "day")
  simdata3$day = 1:62
  simdata3$Date = seq(as.Date("2020-05-01"), by = "day", length.out = 62)
  
  phi = 5.0 # noise level
  dScale = 9.85
  dShape = 1.73
  
  simdata1$Diffs[1] = round(rnbinom(1, size=phi, mu=getlambd(out1,opt_pars1,day=1,ratio = 1.94,
   data = simdata1,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
  simdata1$Cases[1] = 4490 + simdata1$Diffs[1]
  for (i in 2:length(simdata1$day)){
    # gives the model number of observed cases on a given day e.g. our simulated cases per day
    obs = round(rnbinom(1, size=phi, mu=getlambd(out1,opt_pars1,day=i,ratio = 1.94,
                                                 data = simdata1,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
    # we round to the nearest whole number of cases
    simdata1$Cases[i] = simdata1$Cases[i-1] + obs
    simdata1$Diffs[i] = round(obs)
  }
 simdata1
 
 simdata2$Diffs[1] = round(rnbinom(1, size=phi, mu=getlambd(out2,opt_pars2,day=1,ratio = 1.94,
                                                           data = simdata2,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
 simdata2$Cases[1] = 4490 + simdata2$Diffs[1]
 for (i in 2:length(simdata2$day)){
   # gives the model number of observed cases on a given day e.g. our simulated cases per day
   obs = round(rnbinom(1, size=phi, mu=getlambd(out2,opt_pars2,day=i,ratio = 1.94,
                                                data = simdata2,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
   # we round to the nearest whole number of cases
   simdata2$Cases[i] = simdata2$Cases[i-1] + obs
   simdata2$Diffs[i] = round(obs)
 }
 simdata2
 
 simdata3$Diffs[1] = round(rnbinom(1, size=phi, mu=getlambd(out3,opt_pars3,day=1,ratio = 1.94,
                                                            data = simdata2,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
 simdata3$Cases[1] = 4490 + simdata3$Diffs[1]
 for (i in 2:length(simdata3$day)){
   # gives the model number of observed cases on a given day e.g. our simulated cases per day
   obs = round(rnbinom(1, size=phi, mu=getlambd(out3,opt_pars3,day=i,ratio = 1.94,
                                                data = simdata3,sampFrac = 0.35, delayShape = dShape, delayScale = dScale, change_days = "2020-03-14"))  )
   # we round to the nearest whole number of cases
   simdata3$Cases[i] = simdata3$Cases[i-1] + obs
   simdata3$Diffs[i] = round(obs)
 }
 simdata3
 
 #-----------------------------------------------------------
 
 # Non-CI version for testing
 
 times = c(seq(-20.0, max(simdata2$day), 0.1))
 # We use R0 from our pre-SD fits - since we assume this 'baseline' R0 does not change over time. 
 pars = list(N = N, D = 5, R0 = 2.57, k1 = 0.2 ,k2 = 1, q = 0.05, r = 1, ur = 0.4,f = 0.5)
 
 # find the MLE using the data up to and including each day
 mles = c(rep(0.36,dim(simdata2)[1]))
 for (i in 17:dim(simdata3)[1]){
   ffit=optimize(function(f) negloglikef(f, ratio=1.94, pars, simdata2[1:i,], state, times[1:which(times==i)], mysdtiming, overdisp = 5.0), lower=0.0, upper=1.0)
   mles[i] = ffit$minimum
 }
 
 # We want to 'accept' the MLE when it changes by less than 5% over 3 successive days
 for (i in 17:dim(simdata2)[1]){
   if (mles[i]<0.999 && mles[i]>0.0001){ # to stop it accepting right away 
     if (abs((mles[i+1]-mles[i])/mles[i])<0.05){ # compare ith to the i+1th
       if (abs((mles[i+2]-mles[i+1])/mles[i+1])<0.05){
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
 
 # plot the results
 df <- data.frame(Date=simdata2[,1],
                  fhat=mles)
 head(df)
 p<-ggplot(df, aes(x=Date, y=fhat)) + geom_line() + geom_point()
 p <- p + theme_minimal() + scale_color_brewer(palette="Dark2") +
   labs(y = "MLE of f")
 p <- p + geom_vline(xintercept = as.numeric(df$Date[17]), 
                     color = "gray", size=1.5, alpha = 0.5)
 p <- p + annotate(geom = "text", x = df$Date[17]-1, y = min(df$fhat)+0.3, label = "Distancing relaxed", color = "gray", angle = 90) + 
        geom_rect(aes(xmin = df$Date[rec_start], xmax = df$Date[rec_end]), colour = "gray",
                  ymin = 0, ymax = 1, alpha = 0.02) +
   annotate(geom = "text", x = df$Date[rec_start-1], y = 0.75, label = "MLE of f accepted",
            color = "gray", angle = 90) 
 p
 # the lh is completely flat until the 17th, because we don't even consider f until that time (it's
 # not that it thinks it's 1) - but we are assuming it's 0.36 in the model so we set it to that. 
 
 
 #-----------------------------------------------------------
 
 # Now run the daily-MLE method
 
 Nreps = 50
 
 R0_vec = rnorm(Nreps, mean=2.57, sd=0.15)
 
 times = c(seq(-20.0, max(simdata1$day), 0.1))
 
 fit_f = function(state,
                  pars,
                  R0,
                  data1, data2, data3, 
                  ratio){
   library(deSolve)
   pars$R0 = R0
   MLEs1 = rep(1, dim(data1)[1])
   MLEs2 = rep(1, dim(data2)[1])
   MLEs3 = rep(1, dim(data3)[1])
   for( j in 14:dim(data1)[1] ){
     ffit1 = optimize(function(f){ negloglikef(f,
                                              ratio,
                                              pars,
                                              data1[1:j,],
                                              state,
                                              times[1:which(times==j)], mysdtiming, overdisp = 5.0) }, lower=0, upper=1)
     ffit2 = optimize(function(f){ negloglikef(f,
                                               ratio,
                                               pars,
                                               data2[1:j,],
                                               state,
                                               times[1:which(times==j)], mysdtiming, overdisp = 5.0) }, lower=0, upper=1)
     ffit3 = optimize(function(f){ negloglikef(f,
                                               ratio,
                                               pars,
                                               data3[1:j,],
                                               state,
                                               times[1:which(times==j)], mysdtiming, overdisp = 5.0) }, lower=0, upper=1)
     
     MLEs1[j] = ffit1$minimum
     MLEs2[j] = ffit2$minimum
     MLEs3[j] = ffit3$minimum
   }
   
   return( cbind(MLEs1, MLEs2, MLEs3) )
 }
 
 ncore = 4
 cl = makeCluster( ncore )
 clusterExport(cl, varlist=ls(), envir=environment())
 
 system.time({
   results = parLapply( cl, X=1:Nreps, function(i){fit_f(state,
                                                         opt_pars,
                                                         R0_vec[i],
                                                         data1 = simdata1, data2 = simdata2, data3 = simdata3, 
                                                         1.94)})})
 stopCluster(cl)
 
 t = dim(simdata1)[1]
 
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
 for (i in 14:dim(simdata1)[1]){
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
 
 
 save.image("Figure4_relax.RData")


tib2 = res2 %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975))
tib2[1:16,2:4] = c(0.36,0.36,0.36) # during the pre-relaxation period, we are sure that f=0.36.

tib = res %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975)) 
tib[1:16,2:4] = c(0.36,0.36,0.36)

tib3 = res3 %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975)) 
tib3[1:16,2:4] = c(0.36,0.36,0.36)

p1 <- ggplot(data = tib, aes(x=t,y=median))+
   geom_ribbon(data=tib2[1:rec_start2,], aes(x=tib2$t[1:rec_start2], ymin=tib2$L[1:rec_start2], ymax=tib2$U[1:rec_start2]), alpha=0.5, fill="gray50") +
  
  geom_ribbon(data=tib2[-(1:rec_start2-1),], aes(x=tib2$t[-(1:rec_start2-1)], ymin=tib2$L[-(1:rec_start2-1)], ymax=tib2$U[-(1:rec_start2-1)]), alpha=0.5, fill="skyblue1") +
  
  geom_point(data=tib2[1:rec_start2,], x=tib2$t[1:rec_start2],y=tib2$median[1:rec_start2],size=1.6,shape=17, color="gray39") +
  geom_point(data=tib2[-(1:rec_start2-1),], x=tib2$t[-(1:rec_start2-1)],y=tib2$median[-(1:rec_start2-1)],size=1.6,shape=17, color="royalblue3") +
  
  labs(x='Time (days)', y="f (MLE)") +
  theme_minimal() +
  geom_vline(xintercept = as.numeric(res$t[17]), color = "gray25", size=1.5, alpha = 0.5) +
  annotate(geom = "text",
           x = res$t[16]-0.75,
           y = min(res$mle)+0.7,
           label = "Distancing relaxed",
           color = "black", angle = 90) +
  
  annotate("rect", xmin=as.numeric(res2$t[rec_start2]), xmax = as.numeric(res2$t[rec_end2]), 
           ymin=0,ymax=1 ,fill="royalblue3", alpha = 0.5) +
  
  annotate(geom = "text", x = res2$t[rec_start2-1] - 0.75, y = 0.25,
           label = "MLE accepted", color = "royalblue3", angle = 90) 

case_data <- simdata1 %>%
  mutate(Time = as.numeric(Date)-min(as.numeric(Date)))
case_data <- case_data[,c(5,3)]
case_data <- cbind(case_data, simdata2[,3], simdata3[,3])
names(case_data) = c("Time", "Diffs1", "Diffs2", "Diffs3")

p2 <- ggplot(data = case_data, aes(x=Time,y=Diffs1))+
  geom_point(aes(x=Time,y=Diffs2), color = "black")+geom_line(aes(x=Time,y=Diffs2), color = "black") +
  
  labs(x='', y="Daily cases") +theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(breaks=c(0,5,10))
#ggsave('Figure4_casecounts.pdf', width=6, height=3, units = 'in', dpi = 900)

grid.arrange(grobs = list(p2, p1), nrow = 2, heights = c(1,3))
# saved as 6x4 inches


 