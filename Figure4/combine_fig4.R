# Combine figures 4a and 4b to a single pdf
rm(list=ls())
library(ggplot2)
source("../sd_funcs.R")
library(tidyverse)
library(lubridate)
library(gridExtra)

# 4a
load('Figure4.RData')

plot_data = res %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975))


p1_4a <- ggplot(data = plot_data, aes(x=t,y=median))+
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
  theme_minimal() + theme(plot.margin = unit(c(0.2,0,1,0), "cm"))
                    # (add margin for grid arrange)
res
#ggsave('Figure4.pdf', width=6, height=3, units = 'in', dpi = 900)

p2_4a <- bcdata %>%
  mutate(Time = as.numeric(Date)-min(as.numeric(Date))) %>%
  ggplot(aes(x=Time,y=Diffs))+geom_point()+geom_line() + 
  labs(x='', y="Daily cases") +theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
   theme(plot.margin = unit(c(0,0,0,0), "cm"))
# (add margin for grid arrange)


library(gridExtra )
p_4a <- grid.arrange(grobs = list(p2_4a, p1_4a), nrow = 2, heights = c(1,3))

p_4a


#4b
load('Figure4_relax_updated.RData')

tib2 = res %>%
  group_by(t) %>%
  summarise( L=quantile(mle, probs=0.025),
             median=quantile(mle,probs=0.5),
             U = quantile(mle,probs=0.975))
tib2[1:16,2:4] = c(0.36,0.36,0.36) # during the pre-relaxation period, we are 'sure' that f=0.36.


p3 <- ggplot(data = tib2, aes(x=t,y=median))+
  geom_ribbon(data=tib2[1:rec_start2,], aes(x=tib2$t[1:rec_start2], ymin=tib2$L[1:rec_start2], ymax=tib2$U[1:rec_start2]), alpha=0.5, fill="gray50") +
  
  geom_ribbon(data=tib2[-(1:rec_start2-1),], aes(x=tib2$t[-(1:rec_start2-1)], ymin=tib2$L[-(1:rec_start2-1)], ymax=tib2$U[-(1:rec_start2-1)]), alpha=0.5, fill="skyblue1") +
  
  geom_point(data=tib2[1:rec_start2,], x=tib2$t[1:rec_start2],y=tib2$median[1:rec_start2],size=1.6,shape=17, color="gray39") +
  geom_point(data=tib2[-(1:rec_start2-1),], x=tib2$t[-(1:rec_start2-1)],y=tib2$median[-(1:rec_start2-1)],size=1.6,shape=17, color="royalblue3") +
  
  labs(x='Time (days)', y="f (MLE)") +
  theme_minimal() +
  geom_vline(xintercept = as.numeric(res$t[17]), color = "gray25", size=1.5, alpha = 0.5) +
  annotate(geom = "text",
           x = res$t[16]-0.75,
           y = min(res$mle)+0.68,
           label = "Distancing relaxed",
           color = "black", angle = 90) +
  
  annotate("rect", xmin=as.numeric(res$t[rec_start2]), xmax = as.numeric(res$t[rec_end2]), 
           ymin=0,ymax=1 ,fill="royalblue3", alpha = 0.5) +
  
  annotate(geom = "text", x = res$t[rec_start2-1] - 0.75, y = 0.25,
           label = "MLE accepted", color = "royalblue3", angle = 90) + 
  theme(plot.margin = unit(c(0.2,0,1,0), "cm"))
      # (add margin for grid arrange)

case_data <- bcdata %>%
  mutate(Time = as.numeric(Date)-min(as.numeric(Date)))
case_data <- case_data[,c(5,2)]
names(case_data) = c("Time", "Diffs1")

p4 <- ggplot(data = case_data, aes(x=Time,y=Diffs1))+
  geom_point(aes(x=Time,y=Diffs1), color = "black") +
  geom_line(aes(x=Time,y=Diffs1), color = "black") +
  labs(x='', y="Daily cases") +theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(breaks=c(0,10,20,30)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
# (add margin for grid arrange)
#ggsave('Figure4b_multipanel_updated.pdf', width=6, height=3, units = 'in', dpi = 900)

p_4b <- grid.arrange(grobs = list(p4, p3), nrow = 2, heights = c(1,3))
# saved as 6x4 inches

p_4b



# Combine to one figure, margin space already added
grid.arrange(grobs = list(p2_4a, p1_4a, p4, p3), nrow = 4, heights = c(1,3, 1, 3))

# printed 6x4.8inches
