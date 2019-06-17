
# Estimate multivariate linear model of data.
#
# @param data multivariate dataset (matrix).
# @param settings algorithm settings.
# @return set of linear models.

trainMultivariateModel <- function(data, settings){
  # resulting array of models
  models <- matrix(0.0, nrow = settings$q, ncol = length(data[1,])-settings$q+1)
  #all x from (q+1) column to the end of the matrix
  X <- data[,(settings$q+1):ncol(data)]

  for (i in 1:settings$q)
  {
    # i-th column is i-th response
    Y <- data[,i]
    model <- l1fit(X, Y, intercept = TRUE)
    models[i,] <- model$coefficients
  }

  return (models)
}


# Estimate univariate linear model of data.
#
# @param data multivariate dataset (matrix).
# @param i the index of a response to be analyzed.
# @param settings algorithm settings.
# @return a linear model
trainUnivariateModel <- function(data,i,settings){
  if (missing(i)){i=1}
  Y <- data[,i] #i-th column is y
  X <- data[,(settings$q+1):ncol(data)]

  model <- l1fit(X, Y, intercept = TRUE)

  return (model)
}

# Calculate a vector of multivariate model residuals.
#
# @param data multivariate dataset (matrix).
# @param models matrix of linear models coefficients.
# @param settings algorithm settings.
# @return resulting vector of residuals
getResidualsMultivariate <-function(data, models, settings){
  matrixOfResiduals <- getMatrixOfResiduals(data, models, settings)
  result <- vector()
  for (i in 1:length(matrixOfResiduals[1,]))
  {
    result[i] <- sum(abs(matrixOfResiduals[,i])^2)^(1/2)#sum(matrixOfResiduals[,i])
  }
  return (result)
}

# Calculate a vector of univariate model residuals.
#
# @param data multivariate dataset (matrix).
# @param model vactor of linear models coefficients.
# @param settings algorithm settings.
# @param i the index of a response to be analyzed.
# @return resulting vector of residuals
getResidualsUnivariate <-function(data,model,i,settings){
  if (missing(i)){i=1}
  # get independent variables
  X <- data[,(settings$q+1):ncol(data)]

  # multiply weights by variables
  predictions <- as.matrix(X) %*% model$coefficients[2:(length(data[1,])-settings$q+1)]  + model$coefficients[1]

  # calculate absolute residuals
  result <- abs(data[,i]-predictions)

  return(result)

}


# Calculate matrix of multivariate model residuals.
#
# @param data multivariate dataset (matrix).
# @param models matrix of linear models coefficients
# @param settings algorithm settings.
# @return resulting matrix of residuals
getMatrixOfResiduals <-function(data, models, settings){
  # prepare the matrix variable
  result <- matrix(0.0, nrow = settings$q, ncol = length(data[,1]))
  # all x from (q+1) column to the end of the matrix
  X <- data[,(settings$q+1):ncol(data)]

  for (i in 1:settings$q)
  {
    #multiply weights by variables
    predictions <- as.matrix(X) %*% models[i,2:length(models[1,])]  + models[i,1]

    #calculate absolute residuals
    result[i,] <- abs(data[,i]-predictions)
  }

  return (result)
}

# Estimate variable indexes for the given dataset.
#
# @param data multivariate dataset (matrix).
# @param settings algorithm settings.
# @param index index of y if univariate estimation needed (default value is -1)
# @return resulting list of variables
getVariableIndexes <- function(data,settings,index = -1){

  # extract x matrix
  X <- as.matrix(data[,(settings$q+1):ncol(data)])

  # check x dimensionality
  dimx = dim(X)
  if (dimx[2] == 1)
  {
    relevantIndexes <- 1
  }
  else
  {
    relevantIndexes <- 0

    # define the responses to analyze
    rangeOfResponses <- c(1:settings$q)

    # if a particular response selected
    if (index!=-1)
      rangeOfResponses <- c(index)

    # a matrix to store coeffs and their indexes for all the responses
    # (index,value)
    allImpacts <- matrix(0, nrow = settings$l*length(rangeOfResponses), ncol = 2)

    # current position in the matrix
    currentsize <- 0
    # iterate over all the responses
    for (i in rangeOfResponses)
    {
      # extract single response (i.e. i-th column (i-th y))
      Y <- data[,i]
      # run Lasso
      glmnet1 <- cv.glmnet(X,Y)

      # estimated coefficients
      c<-coef(glmnet1,s='lambda.min',exact=TRUE)

      # get rid of bias
      c<-c[-1]

      # get boundary for non-significant contribution within
      # the current response
      empiricalBoundary <- abs(c[which.max(abs(c))])*0.01

      # shrink to zero non-significant coeffs
      # (using max(l, empiricalBoundary))

      # sort absolute contributions of all vars
      orderedArray <- sort(abs(c), decreasing=T)
      restrictedBoundary <- orderedArray[settings$l]

      # final boundary for shrinking
      boundary <- max(abs(empiricalBoundary),abs(restrictedBoundary))

      if (boundary!=0){
        # first round selection
        indexOfMax <- which.max(abs(c))
        relevantIndexes <- c(relevantIndexes,indexOfMax)

        # all relevant indexes
        indexes <- which(abs(c) >= boundary, arr.ind=TRUE)
        values <- c[abs(c)>= boundary]

        # update resulting matrix
        size <- length(indexes)
        allImpacts[(currentsize+1):(currentsize+size),1]<-indexes
        allImpacts[(currentsize+1):(currentsize+size),2]<-values
        currentsize <- size
      }

    }

    # clear definitely relevant indexes
    relevantIndexes<-relevantIndexes[-1]
    relevantIndexes <- unique(relevantIndexes)

    # sort by importance (descreasing), i.e. sort by 2nd column
    allImpacts <- allImpacts[order(abs(allImpacts[,2]), decreasing = TRUE),]

    # find the rest relevant indexes
    for (i in 1:length(allImpacts[,1]))
    {
      # if the current index is not yet selected
      # and it is not empty (zero)
      # and we haven't exceeded the limit l
      if ( !any(relevantIndexes == allImpacts[i,1]) && allImpacts[i,1]!=0 && length(relevantIndexes)<settings$l)
        relevantIndexes <- c(relevantIndexes,allImpacts[i,1])
    }

    # sort final result
    relevantIndexes <- sort(relevantIndexes)
  }

  return (relevantIndexes)
}






