###############################################################################
### Required functions for"When do we see the impact of social distancing?" ###
###############################################################################



#############################################################
#  This is the S E1 E2 I Q R with social distancing ODE model
#
# Summary: - Can infect from E2 (presymptomatic) and I (symptomatic) categories
#          - 2 versions of each category - 1 for distancing individuals and 1 for non-distancing
#          - Individuals move in and out of distancing with rates r and ur, respectively
#          - Distancing individuals have contacts reduced by some fraction f
#          - Some individuals are quarantined and stop being able to infect - at rate q
#          - Average length of symptomatic infectious period D
socdistmodel <- function(t,
                         state,
                         pars,
                         sdtiming){
  with(as.list(c(state,
                 pars)),
       {
         f = sdtiming(t, pars)
         dSdt = -(R0 / (D+1/k2)) * (I + E2 + f * (Id+E2d)) * S/N - r*S + ur*Sd
         dE1dt = (R0 / (D+1/k2)) * (I + E2 + f * (Id+E2d)) * S/N - k1*E1 -r*E1 + ur*E1d
         dE2dt = k1*E1 - k2*E2 - r*E2 + ur*E2d
         dIdt =  k2*E2 - q*I - I/D - r*I + ur*Id
         dQdt = q*I - Q/D - r*Q + ur*Qd
         dRdt = I/D + Q/D - r*R + ur*Rd

         dSddt = -( f*R0/(D+1/k2)) * (I+E2 + f*(Id+E2d)) * Sd/N + r*S - ur*Sd
         dE1ddt = ( f*R0/(D+1/k2)) * (I+E2 + f*(Id+E2d)) * Sd/N - k1*E1d +r*E1 - ur*E1d
         dE2ddt = k1*E1d - k2*E2d + r*E2 - ur*E2d
         dIddt = k2*E2d - q*Id - Id/D + r*I - ur*Id
         dQddt = q*Id - Qd/D + r*Q - ur*Qd
         dRddt = Id/D + Qd/D + r*R - ur*Rd
           list(c(dSdt,
                dE1dt,
                dE2dt,
                dIdt,
                dQdt,
                dRdt,
                dSddt,
                dE1ddt,
                dE2ddt,
                dIddt,
                dQddt,
                dRddt ))
       })
}



#############################################################################
# Now for functions for the likelihood model, and the negative log likelihood:
#
# Log-likelihood - a function of the model parameters, the ODE model simulations, and the data.
# We assume that the number of cases at any time x is negative binomially distributed
sdloglike = function(out,
                        pars,
                        ratio, # this is the ratio of the sampling fraction before/after any changes in testing criteria
                        bcdata, change_days = "2020-03-14", overdisp=5.0,
                        delayShape = 1.73,
                        delayScale = 9.85) {
  # For each of the days in bcdata, this returns the likelihood of the data given the simulations from the model parameters
  ld = vapply(bcdata$day,
              function(x) getlambd(out, # (see below for function getlambd) - the expected number of cases on day x
                                   pars,
                                   x,
                                   ratio,
                                   data=bcdata,
                                   change_days = change_days, delayShape = delayShape,
                                   delayScale = delayScale),
              FUN.VALUE = 1)
  # returns the log likelihood of seeing the data:
  return( sum(dnbinom( bcdata$Diffs, size = overdisp, mu=ld, log=TRUE)) )
}
#
#
# This function returns the model prediction for all the days in the data, for a given set if model parameters
modelpreds = function(out,
                      pars,
                      ratio = 2,
                      data, sampFrac = 0.35, change_days = "2020-03-14") {
  ld = vapply(data$day,
              function(x) getlambd(out,
                                   pars,
                                   x,
                                   ratio, sampFrac = sampFrac, change_days = change_days),
              FUN.VALUE = 1)
  return(data.frame(day = data$day,
                    prediction = ld,
                    data = data$Diffs))
}
#
#
# This is the function works out the model's predicted predicted number of reported cased per day (e.g. it takes into
# account the delay between symptom onset and being reported as a case)
getlambd = function(out,
                    pars,
                    day,
                    ratio = 2,
                    data = bcdata,
                    sampFrac = 0.35,
                    delayShape = 1.73,
                    delayScale = 9.85,
                    change_days = "2020-03-14") {
  # delayShape and delayScale (of a Weibull distribution) were fit for BC in Anderson et al
  # - this is time delay between symptom onset and reporting
  # change_day are the days on which the testing regime changes

  meanDelay = delayScale * gamma(1 + 1 / delayShape)   # comes out to mean = 8.78 with parameters 1.73, 9.85

  # some error messages:
  try(if(var(diff(out$time)) > 0.005) { stop("approx integral assumes equal time steps")})
  try(if(max(out$time) < day) {stop("model simulation is not long enough for the data")})
  try(if(min(out$time) > day - (2*meanDelay+1)) {stop("we need an earlier start time for the model")})

  # relevant times to identify new cases
  ii = which(out$time > day - 2 * meanDelay & out$time <= day)
  dx = out$time[ii[2]] - out$time[ii[1]] # assumes equal time intervals

  # all new cases arising at each of those times
  incoming = with(pars,
                  { k2 * (out$E2[ii] + out$E2d[ii])})

  # take care of changing testing regime over time - this updated version takes a list of days
  # when the testing regime changed ('change_days') and an accompanying list of how the fraction of
  # cases sampled changed on that day ('ratio'), relative to the first days in the data.
  change_day_indicator = sort(c(min(data$day),data$day[match(as.Date(change_days), data$Date)]))
  if (length(change_day_indicator)==1){change_day_indicator = c(0,1)} # case where change day is before the data
  which.sampling = findInterval(day, change_day_indicator)
  ratio = c(1, ratio)
  thisSamp = sampFrac*ratio[which.sampling]

  # each of the past times' contribution to this day's case count
  ft = thisSamp * incoming * dweibull(x = max(out$time[ii]) - out$time[ii],
                                      shape = delayShape,
                                      scale = delayScale)
  # the model's predicted number of reported cases on some day r is comprised of contributions from previous days
  # weighted by the delay

  # return numerical integral of ft (trapezium rule)
  return(thisSamp * 0.5 * (dx) * (ft[1] + 2*sum(ft[2:(length(ft)-1)]) + ft[length(ft)]))
}
#
#
# Lastly, we define functions to return the negative log likelihood - without needing to run the ODEs first
# Firstly, as a function of R0
negloglike = function(R0, ratio, pars, bcdata, state, times, mysdtiming,
                      change_days = "2020-03-14", overdisp = 5.0,
                      delayShape = 1.73, delayScale = 9.85) {

  pars$R0 = R0
  out = as.data.frame(ode(y = state,
                          times = times,
                          func = socdistmodel,
                          parms = pars,
                          sdtiming=mysdtiming)) # This runs the ODE system, starting at 'state' and for parameters 'pars'

  return(-sdloglike(out, pars,ratio,bcdata, change_days = change_days, overdisp = overdisp,
                    delayShape = delayShape, delayScale = delayScale))
}
#
# Then, as a function of f (strength of SD)
negloglikef = function(f, ratio, ode_pars, bcdata, state, times, mysdtiming,
                       change_days = "2020-03-14", overdisp = 5.0,
                       delayShape = 1.73, delayScale = 9.85) {

  ode_pars$f = f
  out = as.data.frame(ode(y = state,
                          times = times,
                          func = socdistmodel,
                          ode_pars,
                          sdtiming=mysdtiming)) #mysdtiming is the function which determines when distancing is turned on/off

  return(-sdloglike(out,ode_pars,ratio,bcdata, change_days = change_days, overdisp = overdisp,
                    delayShape = delayShape, delayScale = delayScale))
}



#########################################################################################################
# This function runs the ODE model lots of times for bootstrapped R0 (normally distributed about our MLE)
#
#
multisolve = function(params,
                      state, times,
                      nReps,
                      timefunc) {
  # get nReps random samples of R0
  rs = rnorm(nReps,
             mean = params$R0,
             sd=0.05)
  # For each sample, runs the ODE model
  biglist = lapply(rs,
                   function(x) {
                     thispars=params
                     thispars$R0=x
                     return( as.data.frame(deSolve::ode(y = state,
                                               times = times,
                                               func = socdistmodel,
                                               parms = thispars,
                                               sdtiming = timefunc)))})
  names(biglist) = rs
  return(dplyr::bind_rows(biglist,
                          .id="R0"))
}



##################################################################################
# These functions get us the number of forecasted cases over time, with error bars
#
# First, just those in the I and Id categories
getCasesbyDay2 = function(df,
                          times,
                          nReps,
                          startDate=dmy("12-03-2020")) {
  CaseInfo=t(vapply(times,
                    function(x) {
                      wherenow = which(df$time==x)
                      normals = df$I[wherenow]
                      selfisols = df$Id[wherenow]
                      return(quantile(normals+selfisols,
                                      probs = c(0.25,0.5,0.75)))
                    },
                    FUN.VALUE = c(2,3,4)))
  return(data.frame(times = times,
                    dates = startDate + times,
                    lower25 = CaseInfo[,1],
                    median = CaseInfo[,2],
                    upper75 = CaseInfo[,3]))
}
#
# Then, all infectives and exposeds (I, Id, E1,E1d, E2, E2d)
getAllCasesbyDay2 = function(df,
                             times,
                             nReps,
                             startDate=dmy("12-03-2020")) {
  CaseInfo=t(vapply(times,
                    function(x) {
                      wherenow = which(df$time == x)
                      normals = (df$I + df$E1 + df$E2)[wherenow]
                      selfisols=(df$Id + df$E1d + df$E2d)[wherenow]
                      return(quantile(normals + selfisols,
                                      probs = c(0.25,0.5,0.75)))
                    },
                    FUN.VALUE = c(2,3,4)))
  return(data.frame(times = times,
                    dates = startDate + times,
                    lower25 = CaseInfo[,1],
                    median = CaseInfo[,2],
                    upper75 = CaseInfo[,3]))
}



##########################################################################
# Plotting function: take in 4 solutions in list form, and make the ggplots
#
makePlots = function(tt2,
                     type="symp",
                     onlyafter = 5,
                     PopScale = TRUE,
                     popSize = N, times) {
  if(type == "symp") {
    cbd1 = getCasesbyDay2(tt2[[1]],
                          times,
                          nReps)
    cbd2 = getCasesbyDay2(tt2[[2]],
                          times,
                          nReps)
    cbd3 = getCasesbyDay2(tt2[[3]],
                          times,
                          nReps)
    cbd4 = getCasesbyDay2(tt2[[4]],
                          times,
                          nReps)
  }
  if(type=="all") {
    cbd1 = getAllCasesbyDay2(tt2[[1]],
                             times,
                             nReps)
    cbd2 = getAllCasesbyDay2(tt2[[2]],
                             times,
                             nReps)
    cbd3 = getAllCasesbyDay2(tt2[[3]],
                             times,
                             nReps)
    cbd4 = getAllCasesbyDay2(tt2[[4]],
                             times,
                             nReps)
  }
  if(PopScale ==TRUE) {
    cbd1[, 3:5] = cbd1[, 3:5] / popSize
    cbd2[, 3:5] = cbd2[, 3:5] / popSize
    cbd3[, 3:5] = cbd3[, 3:5] / popSize
    cbd4[, 3:5] = cbd4[, 3:5] / popSize
  }

  p1 = ggplot(data=cbd1) +
    geom_line(aes(x = dates,
                  y = median)) +
    geom_ribbon(aes(x = dates,
                    ymin = lower25,
                    ymax = upper75),
                alpha = 0.5,
                fill="grey") +
    geom_line(data = cbd2,
              aes(x = dates,
                  y = median)) +
    geom_ribbon(data = cbd2,
                aes(x = dates,
                    ymin = lower25,
                    ymax = upper75),
                alpha = 0.3,
                fill="blue") +
    geom_line(data=cbd3,
              aes(x = dates,
                  y = median)) +
    geom_ribbon(data = cbd3,
                aes(x = dates,
                    ymin = lower25, ymax = upper75),
                alpha = 0.3,
                fill="green") +
    geom_line(data=cbd4,
              aes(x = dates,
                  y = median)) +
    geom_ribbon(data = cbd4,
                aes(x = dates,
                    ymin = lower25, ymax = upper75),
                alpha = 0.3,
                fill="orange")
  p2 = ggplot(data=dplyr::filter(cbd1,
                                 times >= onlyafter)) +
    geom_line(aes(x = dates,
                  y = median))+
    geom_ribbon(aes(x = dates,
                    ymin = lower25,
                    ymax = upper75),
                alpha = 0.5,
                fill = "grey") +
    geom_line(data = dplyr::filter(cbd2,
                                   times >= onlyafter),
              aes(x = dates,
                  y=median)) +
    geom_ribbon(data = dplyr::filter(cbd2,
                                     times >= onlyafter),
                aes(x = dates,
                    ymin = lower25,
                    ymax = upper75),
                alpha = 0.3, fill="blue") +
    geom_line(data = dplyr::filter(cbd3,
                                   times >= onlyafter),
              aes(x = dates,
                  y = median)) +
    geom_ribbon(data = dplyr::filter(cbd3,
                                     times >= onlyafter),
                aes(x = dates,
                    ymin = lower25,
                    ymax = upper75),
                alpha = 0.3,
                fill="green") +
    geom_line(data = dplyr::filter(cbd4,
                                   times >= onlyafter),
              aes(x = dates,
                  y = median)) +
    geom_ribbon(data = dplyr::filter(cbd4,
                                     times >= onlyafter),
                aes(x = dates,
                    ymin = lower25,
                    ymax = upper75),
                alpha = 0.3,
                fill="orange")
  return(list(p1,
              p2))
}



###################################################
# firstdiverge2 by Lloyd T. Elliott
# Precondition: tt a list in the same format as the
#  return value of multisolve. cc is the sensitivity
#  of the test - two time series will be considered
#  different if one exceeds the other by cc cases.
#  The default value of cc is 10.
#
#  The first element of the list tt should be a
#    baseline against which the other elements of tt
#    will be tested.
#
#  If t1 and t2 are provided, only time points in the
#    range t1 <= t <= t2 are considered as possible
#    change points.
#
# Postcondition:
#   Returns a list with named elements 'results' and
#     'posteriors'.
#
#   The value of 'results' is a dataframe with a
#     column 'time' indicating the first time point for
#     which the difference between I+Id values from
#     tt[[1]] exceeds tt[[m]] by cc cases, with posterior
#     probability > 0.95.
#     The column 'posterior' indicates the posterior
#     probability that tt[[1]] exceeds tt[[m]] by cc cases.
#     at the time indicated in 'results' (i.e., 'posterior'
#     must be > 0.95 if it is not missing).
#
#     So, if tt has 4 elements the rows of 'results'
#     have the following names:
#
#       'tt[[1]] -vs- tt[[2]]'
#       'tt[[1]] -vs- tt[[3]]'
#       'tt[[1]] -vs- tt[[4]]'
#
#   If no time point has posterior distribution > 0.95 for
#     a given pair of models, an NA is recorded in
#     'posterior' and also in 'time' columns.
#
#   The value of 'posteriors' indicates the posterior
#     probability that tt[[1]] exceeds tt[[m]] by 10 cases,
#     for EACH time point. (i.e., the first time for which
#     'posteriors' exceeds 0.95 is the time reported in
#     results$time.

#   In all cases, the units of the time points to real
#     times and dates are the same mappings that
#     would be returned from getCasesbyDay2

firstdiverge2 = function(tt, t1, t2, cc = 10, increasing=F) {

  if (cc < 0) {
    warning('firstdiverge: cc must be positive')
    return(NULL)
  }

  # Argument validation

  if (missing(t2)) {
    t2 = max(tt[[1]]$time)
  }

  if (missing(t1)) {
    t1 = min(tt[[1]]$time)
  }

  if (t2 <= t1) {
    warning('firstdiverge: t2 must be after t1')
    return(NULL)
  }

  m = length(tt)
  if (m <= 1) {
    warning('firstdiverge: tt should have more than 1 element.')
    return(NULL)
  }

  support = sort(unique(tt[[1]]$time))
  support = support[support >= t1 & support <= t2]
  n = length(support)

  if (n <= 1) {
    warning('firstdiverge: support of timeseries too small')
  }

  results = data.frame(
    time = rep(NA, m - 1),
    posterior = rep(NA, m - 1)
  )

  for (m1 in 2:m) {
    for (m2 in 1:(m1-1)) {
      if (!all(sort(unique(tt[[m1]]$time)) == sort(unique(tt[[m2]]$time)))) {
        warning('firstdiverge: time alignments mismatch')
        return(NULL)
      }
    }
  }


  posteriors = matrix(NA, n, m - 1)
  colnames(posteriors) = 1:(m-1)
  rownames(posteriors) = support

  for (ii in 2:m) {
    colnames(posteriors)[ii - 1] = sprintf('tt[[1]] -vs- tt[[%d]]', ii)
    rownames(results)[ii - 1] = sprintf('tt[[1]] -vs- tt[[%d]]', ii)
  }

  for (ii in 2:m) {
    ps = c()
    first = -1
    for (t in support) {
      indices = tt[[1]]$time == t
      xx = tt[[1]][indices, 'I'] + tt[[1]][indices, 'Id']

      indices = tt[[ii]]$time == t
      yy = tt[[ii]][indices, 'I'] + tt[[ii]][indices, 'Id']

      numer = 0.0
      denom = 0.0

      for (x in xx) {
        for (y in yy) {
          if ( (increasing && y-x>cc)||(!increasing && x-y>cc) ) {
            numer = numer + 1
          }
          denom = denom + 1
        }
      }

      # For small t, range may be fixed

      if (length(unique(xx)) == 1 || length(unique(yy)) == 1) {
        pv = 0
      } else {
        pv = numer / denom
      }

      if (first == -1 && pv > 0.95) {
        first = t
        results$time[ii - 1] = first
        results$posterior[ii - 1] = pv
      }
      ps = c(ps, pv)
    }
    posteriors[, ii - 1] = ps
  }

  return(list(
    results = results,
    posteriors = posteriors
  ))
}


