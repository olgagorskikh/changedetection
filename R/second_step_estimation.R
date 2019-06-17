recursiveDetection <- function(dataset, s, e, tau, responseIndex,settings){

  # get the center of max energy distance
  maxCenter <- getMaximumCenter(dataset,s,e,tau,responseIndex,settings)

  s <- s + (maxCenter-1)*tau
  e <- s + 2*tau
  #drill down
  updatedtau <- floor(tau*settings$gamma)
  if (updatedtau > 10*settings$l){
    #if there is still enough variables to estimate models well
    #rule: 10 observations per variable (Frank Harrell's book, Regression Modeling Strategies)
    recursiveDetection(dataset,s,e,updatedtau,responseIndex,settings)
  } else {
    #exact detection
    result <- getPreciseChangeLocation(dataset,s,tau,responseIndex,settings)
  }
}

getMaximumCenter<-function(dataset,s,e,tau,responseIndex,settings){
  # data length
  T <- e-s
  # number of periods to observe
  settings$numberOfPeriods <- round(T/tau)
  cat("Number of periods = ",settings$numberOfPeriods,"\n")
  stats = rep(0.0,settings$numberOfPeriods-1)

  # get updated left hand variable indexes
  left = dataset[(1+s):(tau+s),]

  # update variable indexes - only based on left vars
  indexes <- getVariableIndexes(left,settings,responseIndex)
  cat(c("variables: ", paste(indexes,sep=" "),"\n"))

  variables <- c(1:settings$q,indexes+settings$q)


  for (i in 1:(settings$numberOfPeriods-1))
  {
    left <- dataset[((i-1)*tau+1+s):(i*tau+s),variables]
    right <- dataset[(i*tau+1+s):((i+1)*tau+s),variables]

    stats[i] <- calculateEnergyDistance(left, right,settings, responseIndex)

    #cat(c("Periods [",(i-1)*tau+1+s,i*tau+s,"] and [",i*tau+1+s,(i+1)*tau+s,"], StatValue=",stats[i],"\n"))
  }

  maximumIndex <- which.max(stats)
  return (maximumIndex)

}

#point - location of a point indicating small area of a change
getPreciseChangeLocation <-function(dataset,s,tau,responseIndex,settings)
{
  errors = rep(0,2*tau)
  leftBoundary <- s

  cat(c("Left boundary=",leftBoundary,"\n"))

  #update variable indexes
  if (leftBoundary-tau<0) start=1
  else start=leftBoundary-tau
  left = dataset[start:leftBoundary,]

  indexes <- getVariableIndexes(left,settings,responseIndex)
  cat(c("variables: ", paste(indexes,sep=" "),"\n"))

  variables <- c(1:settings$q,indexes+settings$q)

  for (i in leftBoundary:(leftBoundary+2*tau))
  {
    if (i-tau<0) start=1
    else start=i-tau

    left = dataset[start:i,variables]
    right = dataset[(i+1):(3*tau+leftBoundary),variables]

    errors[i-leftBoundary+1] = calculateEnergyDistance(left, right,settings, responseIndex)


  }

  #find the minimum energy distance
  result <- which.max(errors)
  cat(c("Found optimal index within the period ",result,"\n"))
  result = result + s -1
  return (result)
}
