
#' Changeable linear structure of the data with given change points
#'
#' @description The function estimates a set of linear models within the given dataset splitted by the given change points. The models are calculated as L1 regression based on a set of valuable predictors selected by lasso estimator.
#'
#' @param x matrix of regressors with variables in columns and observations in rows
#' @param y matrix of responses with variables in columns and observations in rows
#' @param changes a set of structural change points
#' @param l approximate number of contributing variables (default : overall number of regressors)
#' @return a set of linear models along with corresponding contributing variables indexes
#' @examples
#' T<-60
#' change<-35
#' x<-rnorm(n=T, m=0, sd=1)
#' e<-scale(rt(n=T,3), scale=FALSE)
#' y1<-5*x[1:(change-1)]+e[1:(change-1)]
#' y2<--2*x[change:T]+e[change:T]
#' y<-c(y1,y2)
#'
#' model <- changingModel(x=as.data.frame(x),
#'                        y=as.data.frame(y),
#'                        c(change))
#'
#' @export

changingModel <-function(x,y,changes,l=NULL)
{
  if(is.null(x)||length(x)==0)
    stop("x is missing or empty")
  if(is.null(y)||length(y)==0)
    stop("y is missing or empty")

  #check compatibilities of x and y
  if (nrow(x)!=nrow(y))
    stop("x and y have different number of rows")

  if(is.null(changes)||length(changes)==0)
    stop("array of changes cannot be empty")

  dataset <- cbind(y,x)
  settings <- getDefaultSettings(dataset, q=ncol(y), l=l)
  changes <- c(1,changes,length(dataset[,1]))
  result <- list()
  result$models <- vector()
  result$variables <- vector()

  for (i in 1:(length(changes)-1))
  {
    # current data piece
    datapiece <- dataset[(changes[i]):(changes[i+1]-1),]

    #cat("Period ",i,": [",changes[i],",",changes[i+1]-1,"]\n")

    # get contributing variables
    indexes<-getVariableIndexes(datapiece,settings)
    result$variables[i] <- list(indexes)
    variables <- c(1:settings$q,indexes+settings$q)
    #cat(c("variable indexes: \n", paste(indexes,sep=" "),"\n"))

    # estimate the model
    model <- trainMultivariateModel(datapiece[,variables], settings)
      #trainUnivariateModel(datapiece[,variables],1,settings)
    result$models[i] <- list(model)
    #cat(c("variable weights (start with bias): \n", paste(coef(model),sep=" "),"\n"))

  }

  return (result)

}

#' Moving energy distance
#'
#' @description Estimates energy distance \insertCite{rizzo-szekely10}{changedetection} for each point starting from \code{tau+1} to \code{T-tau}, where \code{T} is a data length. In these terms, energy distance for a point means energy distance between the dataset containing \code{tau} observations to the left and the dataset containing \code{tau} observations to the right of the original point. Hence, we are considering a so-called 'moving frame' of length \code{2tau}. The resulting array shows how the energy distance behaves along the period to analyze.
#' @param x matrix of regressors with variables in columns and observations in rows
#' @param y matrix of responses with variables in columns and observations in rows
#' @param l approximate number of contributing variables (Default : overall number of regressors)
#' @param tau length of a splitting period (Default: l*10, which is dictated by The general rule of thumb \insertCite{Harrell}{changedetection})
#' @param alpha parameter for energy distance formula (default: `1`)
#' @return a list of energy distnce values for pairs of adjacent data segments of length tau (moving-frame construction)
#' @examples
#' T<-60
#' change<-35
#' x<-rnorm(n=T, m=0, sd=1)
#' e<-scale(rt(n=T,3), scale=FALSE)
#' y1<-5*x[1:(change-1)]+e[1:(change-1)]
#' y2<--2*x[change:T]+e[change:T]
#' y<-c(y1,y2)
#'
#' med <- movingEnergyDistance(x=as.data.frame(x),
#'                             y=as.data.frame(y))
#' @references
#'     \insertAllCited
#' @export
movingEnergyDistance <- function(x,y,l=NULL,tau=NULL,alpha=1)
{
  if(is.null(x)||length(x)==0)
    stop("x is missing or empty")
  if(is.null(y)||length(y)==0)
    stop("y is missing or empty")

  #check compatibilities of x and y
  if (nrow(x)!=nrow(y))
    stop("x and y have different number of rows")

  dataset <- cbind(y,x)
  settings <- getDefaultSettings(dataset, q=ncol(y), tau = tau, l=l, alpha=alpha)

  result=vector()

  for (i in (settings$tau+1):(length(dataset[,1])-settings$tau))
  {
    # get contributing variables
    datapiece <- dataset[((i-settings$tau):(i-1)),]

    indexes<-getVariableIndexes(datapiece,settings)
    variables <- c(1:settings$q,indexes+settings$q)

    left = dataset[(i-settings$tau):(i-1),variables]
    right = dataset[(i:(i+settings$tau-1)),variables]
    value = calculateEnergyDistance(left,right,settings)

    result= c(result,value)
  }

  #cat(c(paste(result,sep="\n")))
  return (result)
}

