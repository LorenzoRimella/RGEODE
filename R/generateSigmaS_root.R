#' @title Posterior for sigma (generateSigmaS_root)
#'
#' @description Generate Sigma from the full conditional posterior.
#'
#'
#' @details
#' Using the full conditional posterior computed in the paper, we simulate
#' Sigma, which will be used for the computations of variance-covariance matrix
#' of our data.
#'
#'
#' @param A          array \cr
#'                   The sufficient statistic A.
#'
#' @param Z          array \cr
#'                   The sufficient statistic Z.
#'
#' @param U          array \cr
#'                   The simulated values of U during the iteration.
#'
#' @param  N         int \cr
#'                   The number of data.
#'
#' @param  asigma    double \cr
#'                   The shape parameter of the Gamma prior.
#'
#' @param  bsigma    double \cr
#'                   The rate parameter of the Gamma prior.
#'
#' @param  D         int \cr
#'                   The dimension of data.
#'
#' @param  IND       int_array \cr
#'                   The array with the accepted principal axes.
#'
#'
#' @return \code{generateSigmaS_root} returns a value:
#' \item{sigmaS}{
#'                 array \cr
#'                 the simulation from the full conditional posterior of
#'                 Sigma.
#'           }
#'
#' @note We remark that only the activated principal axes are used to
#'       obtain a simulated sample. Hence we set as default value for the
#'       rejected one the last iterations used.
#'
#'
#'
#' @references
#' \itemize{
#'
#'  \item   [1] Y. Wang, A. Canale, D. Dunson.
#'          "Scalable Geometric Density Estimation" (2016).\cr
#'          The implementation of rgeode_root is inspired to the
#'          Matlab implementation of generateSigmaS_root by Ye Wang.
#' }
#'
#'
#' @author L. Rimella, \email{lorenzo.rimella@hotmail.it}
#'
#' @import stats

  generateSigmaS_root<- function(A, Z, U, N, asigma, bsigma, D, IND)
{
  #*************************************************************************
  #***        Author:      L. Rimella <lorenzo.rimella@hotmail.it>       ***
  #***                                                                   ***
  #***        Supervisors: M. Beccuti                                    ***
  #***                     A. Canale                                     ***
  #***                                                                   ***
  #*************************************************************************

  if(nargs() <8) stop("8 inputs are required.")

  #Update sigma square
  SS     = sum(A-(Z[,IND]^2)%*%(1-U[IND]))
  sigmaS = 1/rgamma(1, shape= asigma+D*N/2, rate= (bsigma+SS/2))

  return(sigmaS)

}
