#----------------------------------------------------------
# preliminaries

rm(list=ls())
library(tidyverse)
library(lubridate)
source('../sd_funcs.R')
library(deSolve)
library(parallel)
library(gridExtra)
set.seed(6202857)

if(file.exists('Figure4_relax_updated.RData')){load('Figure4_relax_updated.RData')}
#----------------------------------------------------------
# Read in the files
bcdata = read.csv("../bc_casecounts0107.csv", header=TRUE)
bcdata[,1] = ymd(bcdata[,1])

# Make this data in the same format we've used elsewhere:
names(bcdata) <- c("Date", "Diffs", "Cumul")
bcdata$day = 1:nrow(bcdata) 
head(bcdata)

# Take a look at the data
plot(bcdata$Date, bcdata$Diffs)

#----------------------------------------------------------
# Various parameters

# DE parameters
N = 5.1e6  #population of BC

# for the initial state of the system
pars = list(N = N, D = 5, R0 = 2.57, k1 = 1/5, k2 = 1, q = 0.05, r = 1, ur = 0.4, f = 0.36)
# Baseline R0 = 2.57 from pre-distancing model fit.
# (and psi = 1.94)

fsi = with(pars,
           r/(r+ur))
nsi = 1 - fsi
i0 = 40 
# Needed to find our 'new' initial state for relaxing distancing
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

# social distancing function
mysdtiming=function(t, pars){
  ifelse( (t > startTime & t< stopTime),
          pars$f,
          0.36)
}
startTime = 17 # start of RELAXING distancing, we assume it's 0.36 before then, as found in Fig 4a
stopTime = 500

# 'state' begins at time -20 (May 1st - 21 days = April 10th), so we need to set 
# it up as such

# What fig4a did was run:
mysdtiming_prev=function(t, pars){
  ifelse( (t > startTime_prev & t< stopTime_prev),
          pars$f,
          1)
}
startTime_prev = 18 # start of distancing
stopTime_prev = 1000
times = c(seq(-20.0, 53, 0.1))

out = as.data.frame(ode(y = state,
                        times = times,
                        func = socdistmodel,
                        pars,
                        sdtiming=mysdtiming_prev))

# In this new figure we need start at time -20 = April 10th
times = c(seq(-20.0, max(bcdata$day), 0.1))
# We initialise with the state of the model on April 10th, from Fig 4a i.e. 'out' above
# Under this old time scale, April 10th = day 41 (because day 1 = March 1st)
state = unlist(out[out$time==41,-1])


#-----------------------------------------------------------

# Non-CI version for testing

# find the MLE using the data up to and including each day
mles = c(rep(0.36,dim(bcdata)[1]))
for (i in 17:dim(bcdata)[1]){
  ffit=optimize(function(f) negloglikef(f, ratio=1.94, pars, bcdata[1:i,], state, times[1:which(times==i)], mysdtiming, overdisp = 5.0), lower=0.0, upper=1.0)
  mles[i] = ffit$minimum
}

# profile likelihood for day 19, say:
#i=19
#x <- seq(0,1, length.out=100)
#for (j in 1:length(x)){
#  y[j] <- negloglikef(x[j], ratio=1.94, pars, bcdata[1:i,], state, times[1:which(times==i)], mysdtiming, overdisp = 5.0)
#}
#plot(x,y)
# Very flat over the first days after distacing, as with 4a.

# We want to 'accept' the MLE when it changes by less than 5% over 3 successive days
for (i in 17:dim(bcdata)[1]){
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
df <- data.frame(Date=bcdata[,1],
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
# 

#-----------------------------------------------------------

# Now run the daily-MLE method

Nreps = 50

# sd = 0.15
R0_vec = rnorm(Nreps, mean=2.57, sd=0.2)

fit_f = function(state,
                 pars,
                 R0,
                 data1, 
                 ratio){
  library(deSolve)
  pars$R0 = R0
  MLEs1 = rep(0.36, dim(data1)[1])
  for( j in 14:dim(data1)[1] ){
    ffit1 = optimize(function(f){ negloglikef(f,
                                              ratio,
                                              pars,
                                              data1[1:j,],
                                              state,
                                              times[1:which(times==j)], mysdtiming, overdisp = 5.0) }, lower=0, upper=1)
    MLEs1[j] = ffit1$minimum
  }
  
  return(MLEs1)
}

ncore = 4
cl = makeCluster( ncore )
clusterExport(cl, varlist=ls(), envir=environment())

system.time({
  results = parLapply( cl, X=1:Nreps, function(i){fit_f(state,
                                                        pars,
                                                        R0_vec[i],
                                                        data1 = bcdata, 
                                                        1.94)})})
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
  if (temp$mle[i]<0.999 && temp$mle[i]>0.0001){ # to stop it accepting right away 
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
rec_start2 = i
rec_end2 = i+2


save.image("Figure4_relax_updated.RData")


tib2 = res %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975))
tib2[1:16,2:4] = c(0.36,0.36,0.36) # during the pre-relaxation period, we are 'sure' that f=0.36.


p1 <- ggplot(data = tib2, aes(x=t,y=median))+
  geom_ribbon(data=tib2[1:rec_start2,], aes(x=tib2$t[1:rec_start2], ymin=tib2$L[1:rec_start2], ymax=tib2$U[1:rec_start2]), alpha=0.5, fill="gray50") +
  
  geom_ribbon(data=tib2[-(1:rec_start2-1),], aes(x=tib2$t[-(1:rec_start2-1)], ymin=tib2$L[-(1:rec_start2-1)], ymax=tib2$U[-(1:rec_start2-1)]), alpha=0.5, fill="skyblue1") +
  
  geom_point(data=tib2[1:rec_start2,], x=tib2$t[1:rec_start2],y=tib2$median[1:rec_start2],size=1.6,shape=17, color="gray39") +
  geom_point(data=tib2[-(1:rec_start2-1),], x=tib2$t[-(1:rec_start2-1)],y=tib2$median[-(1:rec_start2-1)],size=1.6,shape=17, color="royalblue3") +
  
  labs(x='Time (days)', y="f (MLE)") +
  theme_minimal() +
  geom_vline(xintercept = as.numeric(res$t[17]), color = "gray25", size=1.5, alpha = 0.5) +
  annotate(geom = "text",
           x = res$t[16]-0.75,
           y = min(res2$mle)+0.68,
           label = "Distancing relaxed",
           color = "black", angle = 90) +
  
  annotate("rect", xmin=as.numeric(res$t[rec_start2]), xmax = as.numeric(res$t[rec_end2]), 
           ymin=0,ymax=1 ,fill="royalblue3", alpha = 0.5) +
  
  annotate(geom = "text", x = res$t[rec_start2-1] - 0.75, y = 0.25,
           label = "MLE accepted", color = "royalblue3", angle = 90) 

case_data <- bcdata %>%
  mutate(Time = as.numeric(Date)-min(as.numeric(Date)))
case_data <- case_data[,c(5,2)]
names(case_data) = c("Time", "Diffs1")

p2 <- ggplot(data = case_data, aes(x=Time,y=Diffs1))+
  geom_point(aes(x=Time,y=Diffs1), color = "black") +
  geom_line(aes(x=Time,y=Diffs1), color = "black") +
  labs(x='', y="Daily cases") +theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(breaks=c(0,10,20,30))
#ggsave('Figure4b_multipanel_updated.pdf', width=6, height=3, units = 'in', dpi = 900)

grid.arrange(grobs = list(p2, p1), nrow = 2, heights = c(1,3))
# saved as 6x4 inches

# How many days to accept?
as.numeric(res$t[rec_end2]) - 17
