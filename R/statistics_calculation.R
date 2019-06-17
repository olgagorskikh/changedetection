#' A test showing whether two datasets have similar linear structure
#'
#' @description The function performs a nonparametric test showing whether two datasets have similar linear structure or not. The test is based on applying energy distance \insertCite{rizzo-szekely10}{changedetection} to residuals estimated for each dataset separately but only by one model (either first or second). It is implemented as a permutation test with \code{R} rounds and corresponding \code{pzero} \insertCite{gorskikh17}{changedetection}.
#'
#' @param x1 matrix of first period regressors with variables in columns and observations in rows
#' @param y1 matrix of first period responses with variables in columns and observations in rows
#' @param x2 matrix of second period regressors with variables in columns and observations in rows
#' @param y2 matrix of second period responses with variables in columns and observations in rows
#' @param l approximate number of contributing variables (default: overall number of regressors)
#' @param R number of bootstrap rounds (default: `1000`)
#' @param pzero trust level for bootstrap (default: `0.05`)
#' @param alpha parameter for energy distance formula (default: `1`)
#' @return \code{TRUE} or \code{FALSE}
#' @examples
#' T<-60
#' change<-35
#' x<-rnorm(n=T, m=0, sd=1)
#' e<-scale(rt(n=T,3), scale=FALSE)
#' y1<-5*x[1:(change-1)]+e[1:(change-1)]
#' y2<--2*x[change:T]+e[change:T]
#' y<-c(y1,y2)
#'
#' testResult <- changeTest(x1=as.data.frame(x[1:T]),
#'                          y1=as.data.frame(y[1:T]),
#'                          x2=as.data.frame(x[31:T]),
#'                          y2=as.data.frame(y[31:T]),
#'                          R=200)
#'
#' @references
#'     \insertAllCited
#' @export
changeTest<-function(x1, y1, x2, y2,l=NULL, R=1000, pzero=0.05, alpha=1)
{
  if(is.null(x1)||length(x1)==0)
    stop("x1 is missing or empty")
  if(is.null(y1)||length(y1)==0)
    stop("y1 is missing or empty")

  if(is.null(x2)||length(x2)==0)
    stop("x2 is missing or empty")
  if(is.null(y2)||length(y2)==0)
    stop("y2 is missing or empty")

  if (ncol(x1)!=ncol(x2))
    stop("number of regressors in x1 and x2 must be the same")

  if (ncol(y1)!=ncol(y2))
    stop("number of responses in y1 and y2 must be the same")

  #check compatibilities of x and y
  if (nrow(x1)!=nrow(y1)||nrow(x2)!=nrow(y2))
    stop("x1 and y1 or x2 and y2 have different number of rows")

  firstDataset <- cbind(y1,x1)
  secondDataset <- cbind(y2,x2)
  settings <- getDefaultSettings(firstDataset, ncol(y1), l=l, R=R, pzero=pzero, alpha=alpha)

  # if variable selection needed
  if (settings$l<ncol(x1))
  {
    firstVariables <- getVariableIndexes(firstDataset,settings)
    secondVariables <- getVariableIndexes(secondDataset,settings)
    regressorsVariables <- sort(unique(c(firstVariables,secondVariables)))
    allVariables <- c(1:settings$q,regressorsVariables+settings$q)
  }
  else
    {
      allVariables <- c(1:(ncol(x1)+ncol(y1)))
    }

  result <- calculateNPTestResult(firstDataset[,allVariables], secondDataset[,allVariables], settings)
  return (result)

}

#' Energy distance between two datasets
#'
#' Estimate energy distance which is a metric that measures the distance between the distributions
#' of random vectors \insertCite{rizzo-szekely10}{changedetection}.
#' The energy distance is zero if an only if the distributions are identical,
#' otherwise it will diverge.
#'
#' @param x1 matrix of first period regressors with variables in columns and observations in rows
#' @param y1 matrix of first period responses with variables in columns and observations in rows
#' @param x2 matrix of second period regressors with variables in columns and observations in rows
#' @param y2 matrix of second period responses with variables in columns and observations in rows
#' @param l approximate number of contributing variables (default: overall number of regressors)
#' @param alpha parameter for energy distance formula (default: `1`)
#' @return energy distance value
#' @examples
#' T<-60
#' change<-35
#' x<-rnorm(n=T, m=0, sd=1)
#' e<-scale(rt(n=T,3), scale=FALSE)
#' y1<-5*x[1:(change-1)]+e[1:(change-1)]
#' y2<--2*x[change:T]+e[change:T]
#' y<-c(y1,y2)
#'
#' ed <- energyDistance(x1=as.data.frame(x[1:30]),
#'                      y1=as.data.frame(y[1:30]),
#'                      x2=as.data.frame(x[31:T]),
#'                      y2=as.data.frame(y[31:T]))
#'
#' @references
#'     \insertAllCited
#' @export
energyDistance<-function(x1, y1, x2, y2,l=NULL,alpha=1)
{
  if(is.null(x1)||length(x1)==0)
    stop("x1 is missing or empty")
  if(is.null(y1)||length(y1)==0)
    stop("y1 is missing or empty")

  if(is.null(x2)||length(x2)==0)
    stop("x2 is missing or empty")
  if(is.null(y2)||length(y2)==0)
    stop("y2 is missing or empty")

  if (ncol(x1)!=ncol(x2))
    stop("number of regressors in x1 and x2 must be the same")

  if (ncol(y1)!=ncol(y2))
    stop("number of responses in y1 and y2 must be the same")

  #check compatibilities of x and y
  if (nrow(x1)!=nrow(y1)||nrow(x2)!=nrow(y2))
    stop("x1 and y1 or x2 and y2 have different number of rows")

  firstDataset <- cbind(y1,x1)
  secondDataset <- cbind(y2,x2)
  settings <- getDefaultSettings(firstDataset, ncol(y1), l=l, alpha=alpha)

  # if variable selection needed
  if (settings$l<ncol(x1))
  {
    firstVariables <- getVariableIndexes(firstDataset,settings)
    secondVariables <- getVariableIndexes(secondDataset,settings)
    regressorsVariables <- sort(unique(c(firstVariables,secondVariables)))
    allVariables <- c(1:settings$q,regressorsVariables+settings$q)
  }
  else
  {
    allVariables <- c(1:(ncol(x1)+ncol(y1)))
  }

  result <- calculateEnergyDistance(firstDataset[,allVariables], secondDataset[,allVariables], settings)
  return (result)

}

# Calculate whether two datasets have the same structure.
#
# @param firstDataset first dataset (data.frame) having responses as first columns and regressors as all the rest
# @param secondDataset second dataset (data.frame) having responses as first columns and regressors as all the rest
# @param settings settings of the algo (contains )
# @return \code{TRUE} or \code{FALSE}

calculateNPTestResult<-function(firstDataset, secondDataset, settings)
{
  if(is.null(firstDataset)||length(firstDataset)==0||is.null(secondDataset)||length(secondDataset)==0)
    stop("datasets to compare cannot be empty")

  # estimate a set of linear models
  models = trainMultivariateModel(firstDataset, settings)

  firstRes = getResidualsMultivariate(firstDataset, models, settings)
  secondRes = getResidualsMultivariate(secondDataset, models, settings)

  differences = getMatrixOfDifferences(firstRes, secondRes,settings)

  result = calculateStatistics(differences, length(firstDataset[,1]), length(secondDataset[,1]), settings)

  return (result)

}

# Calculate energy distance between two datasets.
#
# @param firstDataset first dataset (data.frame) having responses as first columns and regressors as all the rest
# @param secondDataset second dataset (data.frame) having responses as first columns and regressors as all the rest
# @param i index of a response to analyze (leave empty if you want to consider all responses)
# @return energy distance value

calculateEnergyDistance <- function (firstDataset,secondDataset,settings,i=NULL)
{
  if(is.null(firstDataset)||length(firstDataset)==0||is.null(secondDataset)||length(secondDataset)==0)
    stop("datasets to compare cannot be empty")

  #if all responses considered
  if (is.null(i))
  {
    models = trainMultivariateModel(firstDataset, settings)

    firstRes = getResidualsMultivariate(firstDataset, models, settings)
    secondRes = getResidualsMultivariate(secondDataset, models, settings)
  }
  #if just one of a bunch considered
  else
  {
    model = trainUnivariateModel(firstDataset,i,settings)

    firstRes = getResidualsUnivariate(firstDataset,model,i,settings)
    secondRes = getResidualsUnivariate(secondDataset,model,i,settings)
  }

  differences = getMatrixOfDifferences(firstRes, secondRes,settings)
  indexes = getIndexes(length(firstDataset[,1])+length(secondDataset[,1]), length(firstDataset[,1]), TRUE)
  sums = getArrayOfSums(differences, indexes)

  result = calculateStatisticsValue(sums, length(firstDataset[,1]), length(secondDataset[,1]))

  return (result)

}


# Calculate whether two datasets have the same structure.
# This function is for inner usage within the current file.
#
# @param differences matrix of residuals differences based on the two datasets to be analysed.
# @param firstAmount length of the first dataset.
# @param secondAmount length of the second dataset.
# @return \code{TRUE} or \code{FALSE}

calculateStatistics <- function(differences, firstAmount, secondAmount, settings)
{
  n = length(differences[1,])
  flags = rep(0.0, times = settings$R)

  indexes = getIndexes(n, firstAmount, TRUE)
  sums = getArrayOfSums(differences, indexes)
  value = calculateStatisticsValue(sums, firstAmount,secondAmount)
  flags[1] = value

  for (i in 2:settings$R)
  {
    indexes = getIndexes(n, firstAmount, FALSE)

    sums = getArrayOfSums(differences, indexes)

    value = calculateStatisticsValue(sums, firstAmount, secondAmount)
    flags[i] = value
  }

  initialStat = flags[1]

  flags = sort(flags, decreasing = TRUE)

  boundary = ceiling(settings$R*settings$pzero)

  isLeading = (flags[boundary]<initialStat)


  if (isLeading) {return (TRUE)}
  else {return (FALSE)}

}

# Calculate energy distance between two datasets.
# This function is for inner usage within the current file.
#
# @param points1 first dataset.
# @param points2 second dataset.
# @return energy distance value
#
calculateStatisticsValue <- function(sums, firstAmount, secondAmount)
{
  result = 0.0

  a = (2.0 / (firstAmount * secondAmount));
  b = (firstAmount)^(-2)
  c = (secondAmount)^(-2)

  result = a * sums[3] - b * sums[1] - c * sums[2]

  m = (firstAmount * secondAmount) / (firstAmount + secondAmount)
  result = result * m

  return (result)
}

# Calculate 3D array of sums required for energy distance estimation.
#
# @param x matrix of residuals differences.
# @param indexes set of indexes picked for the estimation.
# @return 3D array of sums
#
getArrayOfSums <- function(x, indexes)
{
  result <- rep(0.0, times = 3)

  sum1 <- 0.0 # X-X
  sum2 <- 0.0 # Y-Y
  sum3 <- 0.0 # X-Y

  n <- length(x[1,])

  #Loop through the columns
  for (j in 1:n){
    #Loop through the rows
    for (i in 1:n){
      #first set column
      if (is.element(j,indexes)){
        if (is.element(i,indexes)){
          sum1 <- sum1 + x[i,j]
        } else {
          sum3 <- sum3 + x[i,j]
        }
      } else {
        #second set column
        if (!is.element(i,indexes)){
          sum2 <- sum2 + x[i,j]
        } else {
          sum3 <- sum3 + x[i,j]
        }
      }
    }
  }

  result[1] <- sum1
  result[2] <- sum2
  result[3] <- sum3

  return (result)
}


# Get the indexes to be picked for the bootstrap sampling
#
# @param n overall number of points
# @param m number of points in the first dataset
# @param flag is \code{TRUE} if no permutations needed, \code{FASLE} otherwise
# @return array of indexes
#
getIndexes <- function (n, m, flag)
{
  indexes = rep(0, times = m)

  # permute: select m points from range 1:n
  if (!flag){
    indexes <- sample.int(n, m)#floor(runif(m, min=1, max=n))
  } else { # no permutations
    indexes <- c(1:m)
  }

  return(indexes)
}

# Calculate matrix of residuals differences required for energy distance estimation.
#
# @param x vector of first residuals.
# @param y vector of second residuals.
# @return matrix of differences between \code{x} and \code{y}
#
getMatrixOfDifferences<-function(x, y, settings)
{
  n <- length(x) + length(y)
  firstAmount <- length(x)

  result <- matrix(0.0,nrow = n, ncol = n)

  for (i in 1:n){

    if (i <= firstAmount){ # First horizontal part
      for (j in 1:n){
        if (j <= firstAmount){ # First vertical part
          result[i,j] <- abs(x[i] - x[j])^settings$alpha/2
        } else {
          result[i,j] <- abs(x[i] - y[j-firstAmount])^settings$alpha/2
        }
      }
    } else { # Second horizontal part
      for (j in 1:n){
        if (j <= firstAmount){ # First vertical part
            result[i,j] <- abs(x[j] - y[i - firstAmount])^settings$alpha/2
        } else {
          result[i,j] <- abs(y[j-firstAmount] - y[i-firstAmount])^settings$alpha/2
        }
      }
    }
  }


  return (result)
}
