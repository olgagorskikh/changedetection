#' Sample dataset
#'
#' This is artificial data following changing multivariate liner structure.
#'
#' @format A data frame with 600 rows and 13 columns:
#' \itemize{
#'   \item{\eqn{y_1, y_2, y_3}}{ response variables}
#'   \item{\eqn{x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9, x_10}}{ predictors}
#' }
#'
#'
#' @details
#' The underlying model is divided into 3 regimes:
#' \tabular{llllllll}{
#' 1) regime 1 (rows 1-59):\cr
#'       \tab \eqn{y_1 =  x_1 +  x_5}\cr
#'       \tab \eqn{y_2 =  2x_1 + 7x_9}\cr
#'       \tab \eqn{y_3 = -4x_9}\cr
#' 2) regime 2 (rows 60-279):\cr
#'       \tab \eqn{y_1 =  -3x_2 + x_9}\cr
#'       \tab \eqn{y_2 =  2x_1 + 7x_9}\cr
#'       \tab \eqn{y_3 = -4x_9}\cr
#' 3) for rows 280-400:\cr
#'       \tab \eqn{y_1 =  -3x_2 + x_9}\cr
#'       \tab \eqn{y_2 = 6x_2}\cr
#'       \tab \eqn{y_3 = -4x_9}\cr
#'
#' By construction, the \code{sampledataset} has two structural changes happened at points 60 and 280, which makes it useful to be applied to various functions in the current package.
#'
#' @name sampledataset

"sampledataset"
