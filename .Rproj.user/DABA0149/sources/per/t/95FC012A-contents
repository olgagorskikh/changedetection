# First step estimation. Get the suspected change areas centers
#
# @param dataset A matrix containing both responses and predictors.
# @param settings algorithm settings.
#
# @return The set of approximate change centers indexes \code{flags}.

getInitialFlags <- function(dataset, settings){

  tau <-settings$tau

  flags = rep(FALSE,settings$numberOfPeriods-1)
  stats = rep(0.0,settings$numberOfPeriods-1)

  # assess the statistics flags for all neighbour pairs of periods
  for (i in 1:(settings$numberOfPeriods-1))
  {
    left <- dataset[((i-1)*tau+1):(i*tau),]
    variables<-getVariableIndexes(left,settings)
    cat(c(" variables: ", paste(variables,sep=" "),"\n"))
    if (10*length(variables)>settings$tau)
      stop(c("Number of contributing variables is too large for the chosen tau (minimum length of a regime). You may 1) select tau to be at least ",10*length(variables)," or 2) restrict the number of contributing variables by selecting parameter l to be <= ", settings$tau/10))
    # add responses
    variables<-c(1:settings$q,variables+settings$q)

    left <- dataset[((i-1)*tau+1):(i*tau),variables]
    right <- dataset[(i*tau+1):((i+1)*tau),variables]

    flags[i] <- calculateNPTestResult(left, right, settings)
    stats[i] <- calculateEnergyDistance(left, right, settings)


    #cat(c("Periods [",(i-1)*tau+1,i*tau,"] and [",i*tau+1,(i+1)*tau,"], Flag=",flags[i],"; StatValue=",stats[i],"\n"))
  }

  # exclude excess flags
  for (i in 2:(settings$numberOfPeriods-1)){
    # if a test showed a change in two consequential points
    if(flags[i]&&flags[i-1]){
      # if the current change is stronger than previous
      if (stats[i]> stats[i-1]){
        flags[i-1] <- FALSE
      } else { #otherwise
        flags[i] <- FALSE
      }
    }
  }

  return (flags)
}

# First step estimation. Get the suspected change areas centers
#
# @param dataset A matrix containing both responses and predictors.
# @param changes a set of approximate change points.
# @param settings algorithm settings.
#
# @return The set of response indexes having biggest changes.

getResponseIndexes <- function(dataset,changes,settings)
{
  result <- vector()
  tau <- settings$tau

  for (i in 1:(length(changes)))
  {
    # get updated left hand variable indexes
    left = dataset[(changes[i]*tau - tau):(changes[i]*tau-1),]

    #update valuable variables
    variables <- getVariableIndexes(left,settings)
    variables <- c(1:settings$q,variables+settings$q)

    # cut relevant data pieces
    left = dataset[(changes[i]*tau - tau):(changes[i]*tau-1),variables]
    right = dataset[(changes[i]*tau):(changes[i]*tau+tau),variables]
    stats = rep(0.0,settings$q)

    for (j in 1:settings$q)
    {
      stats[j] = calculateEnergyDistance(left,right,settings,j)
    }

    result <- c(result,which.max(stats))
  }

  return (result)
}

