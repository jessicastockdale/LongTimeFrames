# Combine figures 2a and 2b to a single pdf
rm(list=ls())
library(ggplot2)
source("Figure2bot/sd_funcs.R")
library(tidyverse)
library(lubridate)
library(gridExtra)

# 2b
load('Figure2bot/Figure2.RData')

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

p_2b = last_plot()
p_2b


#2a
load('Figure2top/Figure2.RData')

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

p_2a = last_plot()
p_2a

grid.arrange(grobs = list(p_2a, p_2b), nrow = 2)

# add spacing for a) and b) labels
pl <- list(p_2a, p_2b)
margin = theme(plot.margin = unit(c(0.2,0,1,0), "cm"))
grid.arrange(grobs = lapply(pl, "+", margin))
# printed 6x4.8inches
