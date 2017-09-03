#' @title The core of the Gibbs sampling.
#'
#' @description Perform GEODE without missing data.
#'
#' @param Y          array_like \cr
#'                   a real input matrix (or data frame), with dimensions
#'                   \eqn{(n, D)}. It is the real matrix of data.
#'
#' @param d          int, optional \cr
#'                   it is the conservative upper bound for the dimension
#'                   D. We are confident that the real dimension is smaller
#'                   then it.
#'
#' @param burn       int, optional \cr
#'                   number of burn-in to perform in our Gibbs sampler. It
#'                   represents also the stopping time that stop the choice
#'                   of the principal axes.
#'
#' @param its        int, optional \cr
#'                   number of iterations that must be performed after
#'                   the burn-in.
#'
#' @param tol        double, optional \cr
#'                   threshold for adaptively removing redundant
#'                   dimensions. It is used compared with the ratio:
#'                   \eqn{\frac{\alpha_j^2(t)}{\max \alpha_i^2(t)}}.
#'
#' @param atau       double, optional \cr
#'                   The parameter \eqn{a_\tau} of the truncated
#'                   Exponential (the prior for \eqn{\tau_j}).
#'
#' @param asigma     double, optional \cr
#'                   The shape parameter \eqn{a_\sigma} of the
#'                   truncated Gamma (the prior for \eqn{\sigma^2}).
#'
#' @param bsigma     double, optional \cr
#'                   The rate parameter \eqn{b_\sigma} of the
#'                   truncated Gamma (the prior for \eqn{\sigma^2}).
#'
#'
#' @param starttime  int, optional \cr
#'                   starting time for adaptive pruning. It must be less
#'                   then the number of burn-in.
#'
#'
#' @param stoptime   int, optional \cr
#'                   stop time for adaptive pruning. It must be less
#'                   then the number of burn-in.
#'
#' @param fast       bool, optional \cr
#'                   If \eqn{TRUE} it is run using fast d-rank SVD.
#'                   Otherwise it uses the classical SVD.
#'
#' @param c0         double, optional \cr
#'                   Additive constant for the exponent of the pruning step.
#'
#' @param c1         double, optional \cr
#'                   Multiplicative constant for the exponent of the pruning step.
#'
#' @note The parameters are the same as in \code{rgeode}.
#'
#' @return The outputs are the same as in \code{rgeode}.
#'
#'
#' @references
#' \itemize{
#'
#'  \item   [1] Y. Wang, A. Canale, D. Dunson.
#'          "Scalable Geometric Density Estimation" (2016).\cr
#' }
#'
#'
#' @author L. Rimella, \email{lorenzo.rimella@hotmail.it}
#'
#'
#' @import stats MASS
#'

   rgeode_root<- function(Y, d, burn, its, tol, atau, asigma,
                          bsigma, starttime, stoptime, fast, c0, c1)
{
  #*************************************************************************
  #***        Author:      L. Rimella <lorenzo.rimella@hotmail.it>       ***
  #***                                                                   ***
  #***        Supervisors: M. Beccuti                                    ***
  #***                     A. Canale                                     ***
  #***                                                                   ***
  #*************************************************************************

  # Fit GEODE on dataset with no missing data

  ##########################################################################
  # First step:
  # Computation of the mean and the principal axes and creation
  # of the sufficient statistics A,Z
  ##########################################################################

  #Mean and Principal axes
  N= dim(Y)[1]
  D= dim(Y)[2]

  mu =  matrix(apply(Y,2,mean), D, 1)
  Y_c=  t(apply(Y,1,function(x)return(x-mu)))

  #Perform a fast rank-d SVD or a simple SVD
  if(fast)
  {
    W= randSVD(t(Y_c), k= d )$u
  }

  if(!fast)
  {
    W= svd(t(Y_c), nu= d, nv= d)$u
  }


  #Sufficient Statistics
  A = matrix(apply(Y_c, 1, function(x)return(sum(x^2))), N, 1)
  Z = Y_c %*% W

  ##########################################################################
  # Second step:
  # Preparation for the Gibbs sampler:
  # nc=its nb=burn tol=tol a=atau a_sigma = asigma; b_sigma= bsigma
  # step=step stoptime = ; starttime = ; T= Time
  ##########################################################################

  Time = burn + its

  InD=   rep(list(c(1:d)), stoptime)

  u=      matrix(0.5, d, Time)

  tau=       matrix(1, d, Time)
  tau[,1]= rexptr(d, atau, range= c(1,Inf))

  sigmaS= matrix(1, Time, 1)
  sigmaS[1]= 1/rgamma(1,shape= asigma, rate= bsigma)


  ##########################################################################
  # Third step:
  # Gibbs sampler:
  ##########################################################################

  InDtmp = seq(1,d)

  #summ=0

  for(iter in 2:Time)
  {

    #Update u
    #u[,iter] = generateU_root(Z[,InDtmp], sigmaS[iter-1], N,
    #                          tau[InDtmp,iter-1], u[,iter-1])

    u[,iter]= u[,iter-1]
    u[InDtmp,iter] = CgenerateU_root(tau[InDtmp,iter-1], N, sigmaS[iter-1],
                                     Z[,InDtmp])


    #tic(func.tic = NULL)
    #tau[,iter] = generateTau_root(u[,iter], atau, InDtmp, tau[,iter-1])
    tau[,iter]= tau[,iter-1]
    tau[InDtmp,iter]= CgenerateTau_root(u[InDtmp,iter],tau[InDtmp,iter],
                                        atau, max(InDtmp))
    #c=toc(func.toc = NULL)
    #cc= c$toc- c$tic
    #summ=summ+cc

    sigmaS[iter] = generateSigmaS_root(A, Z, u[,iter], N, asigma, bsigma,
                                       D, InDtmp)

    # Adaptively prune the intrinsic dimension
    # Wang version:
    #if (iter <= stoptime & iter > starttime)
    #{
    #  u_accum = u_accum + ( u[,iter]==1 )
    #  if(any(adptpos == iter))
    #  {
    #    adapt[iter] = 1
    #    ind = InD[[iter-1]]
    #    tmp = u_accum[ind]
    #    vec = (1/u[InD[[iter-1]],iter]-1)*sigmaS[iter]
    #    d1 = ind

    #    if ( sum(tmp > (stoptime-starttime)*tol) + sum(vec/max(vec) < tol) )>0
    #    {
    #      ind = ind[ ind %in% ind[ d1 >= min(d1( (tmp > (stoptime-starttime)*tol) |
    #                                                (vec/max(vec) < tol) )) ] = [];
    #    }
    #    InD[[iter]] = ind

    #    else InD[[iter]]= InD[[iter-1]]
    #

    #  }
    #  InDtmp = InD[[iter]]

    #}

    if(iter <= stoptime & iter > starttime)
    {
      ind= InD[[iter-1]]
      #find our alphas
      vec_alpha= (1/u[ind,iter]-1)*sigmaS[iter]
      max_alpha= max(vec_alpha)
      if(iter== stoptime)
      {
        if( sum(ind[vec_alpha/max_alpha < tol])>0 )
        {
          ind= ind[ ind < min( ind[vec_alpha/max_alpha < tol] ) ]
        }

      }

      else if( runif(1, 0, 1)< p(iter, c0, c1) )
      {
        if( sum(ind[vec_alpha/max_alpha < tol])>0 )
        {
          ind= ind[ ind < min( ind[vec_alpha/max_alpha < tol] ) ]
        }

        else
        {
          if(ind[length(ind)]+1<d)
          {
            ind= c(ind, ind[length(ind)]+1)
          }

        }
      }

      InD[[iter]]= ind

      InDtmp = InD[[iter]]

    }

  }


#print(summ)

return( list(InD, u, tau, sigmaS, W, mu) )
}
