
#' A list of change locations occured in the structure having given responses `y` and regressors `x`
#'
#' @description The main function of the package which estimates a set of structural change points for a dataset following multivariate (or univariate) linear model.
#'
#' We consider the following problem. Assume we have observations \eqn{y_t\in{R}^q}, \eqn{x_t\in{R}^p} over time \eqn{t=1,...,N} which follow linear model
#'
#'                     \deqn{ y_t=  x'_t{\beta}(t)+{\epsilon}_t}
#'
#' where \eqn{\epsilon_t\in{R}^d} stands for the model noise and \eqn{\beta} is a \eqn{p\times d}-dimensional piecewise-constant function, i.e. \eqn{{\beta}(t)\in{R}^{p\times d}} is a weight matrix for every \eqn{t}. A change point is defined as a time instant \eqn{\hat{t}} where \eqn{\beta} shifts. More precisely, a set of points separating sequential regimes is a set of change points of the model and the objective of \code{detectChanges} function is to find this set.
#'
#' The approach is based on the splitting procedure NSA \insertCite{gorskikh17}{changedetection} and energy distance analysis \insertCite{rizzo-szekely10}{changedetection}.
#'
#'
#' @param x matrix of regressors with variables in columns and observations in rows
#' @param y matrix of responses with variables in columns and observations in rows
#' @param l approximate number of contributing variables (default: overall number of regressors)
#' @param tau length of splitting periods (default: l*10, which is dictated by The general rule of thumb \insertCite{Harrell}{changedetection})
#' @param R number of bootstrap rounds (default: `1000`)
#' @param pzero trust level for bootstrap (default: `0.05`)
#' @param gamma tau reduction rate used in NSA (default: `0.5`)
#' @param alpha parameter for energy distance formula (default: `1`)
#' @return A set of change locations indexes \code{changes}.
#' @examples
#' T<-60
#' change<-35
#' x<-rnorm(n=T, m=0, sd=1)
#' e<-scale(rt(n=T,3), scale=FALSE)
#' y1<-5*x[1:(change-1)]+e[1:(change-1)]
#' y2<--2*x[change:T]+e[change:T]
#' y<-c(y1,y2)
#'
#' changes(x=as.data.frame(x),
#'         y=as.data.frame(y),
#'         tau=20, R=100)
#'
#' @references
#'     \insertAllCited
#' @export


changes <-function(x, y, tau=NULL, l=NULL, R=1000, pzero=0.05, gamma=0.5, alpha=1){

  if(is.null(x)||nrow(x)==0)
    stop("x is missing or empty")
  if(is.null(y)||nrow(y)==0)
    stop("y is missing or empty")

  #check compatibilities of x and y
  if (nrow(x)!=nrow(y))
    stop("x and y have different number of rows")

  settings <- getDefaultSettings(cbind(y,x), ncol(y), tau=tau, l=l, R=R, pzero=pzero, gamma=gamma, alpha=alpha)
  changes <- runDetection(cbind(y,x), settings)
  return (changes)

}

# Get a list of change locations, occured in the structure having given responses `y` and regressors `x`
#
# @param dataset a matrix having responses as first columns and regressors as all the rest (the only required argument)
# @param q dimensionality of a response (`1` stands for univariate, `2` for bivariate etc.) (default: `1`)
# @param l approximate number of contributing variables (default : overall number of regressors)
# @param tau length of splitting periods (default: l*10, which is dictated by The general rule of thumb (Frank Harrell, Regression Modeling Strategies))
# @param R number of bootstrap rounds (default: `1000`)
# @param pzero trust level for bootstrap (default: `0.05`)
# @param gamma tau reduction rate (default: `0.5`)
# @param alpha parameter for energy distance (default: `1`)
#
# @return A list containing algo settings \code{settings}.
getDefaultSettings <- function(dataset=NULL, q=NULL, l=NULL, tau=NULL, R=1000, pzero=0.05, gamma=0.5, alpha=1){

  if(is.null(dataset)||length(dataset)==0)
    stop("dataset cannot be empty")

  # Settings initialization
  settings <- list()

  # dimensionality of response
  ifelse(is.null(q), settings$q <- 1, settings$q <- q)

  # number of contributing variables
  ifelse(is.null(l), settings$l <- (ncol(dataset)-settings$q), settings$l <- l)

  # check tau (min length of a period)
  if(is.null(tau)) tau <- 10*settings$l # from The general rule of thumb (Frank Harrell, Regression Modeling Strategies
  if(tau < 1) tau <- floor(nrow(dataset)*tau)
  if(tau <= settings$l)
    stop("minimum period length must be greater than the number of regressors")
  if(tau > floor(nrow(dataset)/2))
    stop("minimum period length must be smaller than half of the number of observations")
  settings$tau <- tau

  # number of bootstrap rounds
  settings$R <- R
  # trust level for bootstrap
  settings$pzero <- pzero
  # tau reduction rate
  settings$gamma <- gamma
  # alpha for energy distance
  settings$alpha <- alpha

  return (settings)
}

runDetection <- function(dataset, settings){

  # data length
  T <- length(dataset[,1])
  # number of periods to observe
  settings$numberOfPeriods <- round(T/settings$tau)
  cat("Number of periods = ",settings$numberOfPeriods,"\n")

  # first step estimation
  flags <- getInitialFlags(dataset,settings)
  settings$numberOfChanges <- length(flags[flags==TRUE])
  cat("Number of changes = ",settings$numberOfChanges,"\n")
  changeIndexes <- which(flags==TRUE)

  if (is.null(length(changeIndexes))||length(changeIndexes)==0){
    cat("There are no changes in the dataset\n")
    return (NULL)
  }

  else{
    cat("Change indexes = ",changeIndexes,"\n")
    responseIndexes <- getResponseIndexes(dataset,changeIndexes,settings)
    cat("Response indexes = ",responseIndexes,"\n")
    # second step estimation
    finalChanges <- rep(0.0,settings$numberOfChanges)

    for (i in 1:settings$numberOfChanges){
      start <- (changeIndexes[i]-1)*settings$tau +1
      end <- (changeIndexes[i]+1)*settings$tau
      cat("Consider period [",start,",",end,"]\n")
      finalChanges[i] <- recursiveDetection(dataset, start, end, floor(settings$tau*settings$gamma), responseIndexes[i],settings)
    }

    cat("Final change points = ",finalChanges,"\n")
    return (finalChanges)
  }
}










