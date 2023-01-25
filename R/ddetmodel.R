#' Density of determinantal point process model via direct calculation
#'
#' Evaluates the density of a determinantal point process model at a given point
#' pattern via direct calculation (i.e. not using FFT, etc.).
#'
#' @param X An object of class [spatstat.geom::ppp] specifying a planar point
#'   pattern.
#' @param model An object of class [spatstat.model::detpointprocfamilyfun]
#'   specifying the model to be used.
#' @param trunc A numeric value specifying how the model truncation is
#'   preformed. Defaults to `0.99`. See details of
#'   [spatstat.model::simulate.detpointprocfamily()].
#' @param log A boolean value specifying whether to compute the logarithm.
#'   Defaults to `TRUE`.
#' @param ... Extra parameters to be passed on to next methods. Ignored.
#'
#' @keywords internal
ddetmodel.direct <- function(X,
                             model,
                             trunc = .99,
                             log = TRUE,
                             ...) {
  if(inherits(model, "dppm"))
    model <- model$fitted
  if(!inherits(model, "detpointprocfamily"))
    stop("Argument model must be of class detmodel.")
  if(!inherits(X, "ppp"))
    stop("Only implemented for two dimensional patterns at the moment. Please
         supply a point pattern of class ppp.")
  W <- spatstat.geom::as.owin(X)
  xyscale <- c(diff(W$xrange), diff(W$yrange))
  XX <- spatstat.geom::affine(X, mat = matrix(c(1/xyscale[1],0,0,1/xyscale[2]), 2, 2),
                              vec = -c(mean(W$xrange), mean(W$yrange))/xyscale)
  n <- XX$n
  x <- XX$x
  y <- XX$y

  ########### Old code in pure R and not exploiting stationarity  ##############
  ## lambda.k <- lambda.and.k(model, trunc, xyscale)
  ## trunc <- lambda.k$trunc
  ## prec <- lambda.k$prec
  ## n <- length(lambda.k$lambda)

  ## if(any(lambda.k$lambda<0|lambda.k$lambda>=1)) return(ifelse(log,-Inf,0))
  ## D <- -sum(log(1-lambda.k$lambda))
  ## lambdat <- lambda.k$lambda/(1-lambda.k$lambda)

  ## basis <- fourierbasis(matrix(c(x,y),ncol=2), lambda.k$k, win = owin(c(-.5,.5),c(-.5,.5))) * sqrt(lambdat)
  ## C <- Re(crossprod(Conj(basis),basis))
  ##############################################################################

  ########### New code using some C and exploiting stationarity  ###############
  lambda.k <- lambda.and.k(model, trunc, xyscale, stationary=TRUE)
  lambdao <- lambda.k$lambdao
  lambdaa <- lambda.k$lambdaa
  lambda <- lambda.k$lambda
  k1a <- lambda.k$k1a
  k2a <- lambda.k$k2a
  trunc <- lambda.k$trunc
  prec <- lambda.k$prec
  rm(lambda.k)

  if(lambdao<0|lambdao>=1|any(lambdaa<0|lambdaa>=1)|any(lambda<0|lambda>=1))
    return(ifelse(log,-Inf,0))
  D <- log(1-lambdao)+2*sum(log(1-lambdaa))+4*sum(log(1-lambda))
  lambdao <- lambdao/(1-lambdao)
  lambdaa <- lambdaa/(1-lambdaa)
  lambda <- lambda/(1-lambda)

  C <- CtildeStat2d_cpp(x,y,lambdao,lambdaa,lambda,k1a,k2a)

  tmp <- determinant(C, logarithm = TRUE)
  logdet <- as.numeric(tmp$modulus)
  detsign <- as.numeric(tmp$sign)

  res <- logdet + D + spatstat.geom::volume(W) - n * log(spatstat.geom::volume(W))

  if(!log)
    res <- exp(res)
  attr(res, "determinantal") <- list(logdet = logdet, D = D, detsign = detsign, trunc = trunc, prec = prec)
  return(res)
  ### numeric with attributes containing intermediate calculations etc.
}

lambda.and.k <- function#Internal function calculating lambda and k
### This function is mainly for internal package use and is usually
### not called by the user.
(model,
 ### object of class \code{"detmodel"}
 trunc,
 ### numeric giving the truncation
 Wscale,
 ### numeric giving the scale of the window relative to a unit box
 stationary = FALSE
 ### logical indicating whether the stationarity of the model should be used (only works in dimension 2).
){
  dim <- spatstat.model::dim.detpointprocfamily(model)
  if(stationary&&dim!=2)
    stop("Stationarity can only be exploited in dimension 2 at the moment.")

  ## Calculate expected number of points if the intensity is a parameter
  expnum <- NULL
  rhoname <- model$intensity
  if(!is.null(rhoname))
    expnum <- getElement(model$fixedpar, rhoname)*prod(Wscale)
  ## Get the maximal truncation in each dimension
  maxtrunc <- 2^14# dppspatstat.options("max.trunc")^(1/dim)
  ## Extract spectral density
  specden <- spatstat.model::dppspecden(model)
  truncrange <- spatstat.model::dppspecdenrange(model) * max(Wscale)

  if(trunc>=1){ ## Integer truncation fixed by user.
    if(stationary){
      ## Coordinates on axes:
      k1a <- c(rep(0,trunc),1:trunc)
      k2a <- c(1:trunc,rep(0,trunc))
      ## Coordinates of ordinary points:
      k1 <- rep(1:trunc,trunc)
      k2 <- rep(1:trunc,each=trunc)
      ## Spectral densities:
      lambdao <- specden(0)
      lambdaa <- specden(sqrt((k1a/Wscale[1])^2+(k2a/Wscale[2])^2))
      lambda <- specden(sqrt((k1/Wscale[1])^2+(k2/Wscale[2])^2))
      prec <- (lambdao+2*sum(lambdaa)+4*sum(lambda))/expnum
    } else{
      trunc <- floor(trunc)
      k <- do.call(expand.grid, replicate(dim, seq(-trunc,trunc), simplify=FALSE))
      kscaled <- k*matrix(1/Wscale, nrow(k), ncol(k), byrow = TRUE)
      if(model$isotropic){
        lambda <- specden(sqrt(rowSums(kscaled^2)))
      } else{
        lambda <- specden(kscaled)
      }
      prec <- sum(lambda)/expnum
    }
  } else{ ## Integer truncation calculated from user-specified precision.
    if(is.null(expnum))
      stop("Cannot calculate truncation adaptively in a model without intensity parameter. Please specify trunc directly as a positive integer.")
    prec0 <- trunc
    trunc <- 1
    prec <- 0
    ## cat("truncation is being calculated adaptively. Current truncation:\n")
    while(prec<=prec0 & (2*trunc)<=maxtrunc & trunc<=truncrange){
      # while(prec<=prec0 & (trunc + 1)<=maxtrunc & trunc<=truncrange){
      trunc <- 2*trunc
      # trunc <- trunc + 1
      if(stationary){
        ## Coordinates on axes:
        k1a <- c(rep(0,trunc),1:trunc)
        k2a <- c(1:trunc,rep(0,trunc))
        ## Coordinates of ordinary points:
        k1 <- rep(1:trunc,trunc)
        k2 <- rep(1:trunc,each=trunc)
        ## Spectral densities:
        lambdao <- specden(0)
        lambdaa <- specden(sqrt((k1a/Wscale[1])^2+(k2a/Wscale[2])^2))
        lambda <- specden(sqrt((k1/Wscale[1])^2+(k2/Wscale[2])^2))
        prec <- (lambdao+2*sum(lambdaa)+4*sum(lambda))/expnum
      } else{
        k <- do.call(expand.grid, replicate(dim, seq(-trunc,trunc), simplify=FALSE))
        kscaled <- k*matrix(1/Wscale, nrow(k), ncol(k), byrow = TRUE)
        if(model$isotropic){
          lambda <- specden(sqrt(rowSums(kscaled^2)))
        } else{
          lambda <- specden(kscaled)
        }
        prec <- sum(lambda)/expnum
      }
    }
    ## cat("\n")
    if(prec<prec0){
      warning(paste0("Adaptive truncation stopped at ", trunc, ". The precision is only ", prec))
    }
  }
  if(stationary){
    rslt <- list(lambdao=lambdao, lambdaa=lambdaa, lambda=lambda, k1a=k1a, k2a=k2a)
  } else{
    rslt <- list(lambda=lambda, k=k)
  }
  return(c(rslt, list(prec=prec, trunc=trunc)))
  ##value<< A list
}
