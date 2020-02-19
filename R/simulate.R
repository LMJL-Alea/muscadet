#' @export
rbidppn <- function(n = 1, seed = 1234, rho1 = 100, rho2 = 100, tau = 0.2,
                    alpha1 = 0.03, alpha2 = 0.03, alpha12 = 0.05,
                    nu1 = 10, nu2 = 10, nu12 = 10,
                    progress = TRUE,
                    Kspec="Kspecmatern",
                    testtau="testtaumatern")
{
  set.seed(seed)
  future::plan(future::multiprocess)
  furrr::future_map(
    .x = 1:n,
    .f = ~ rbidpp(
      rho1 = rho1, rho2 = rho2, tau = tau,
      alpha1 = alpha1, alpha2 = alpha2, alpha12 = alpha12,
      nu1 = nu1, nu2 = nu2, nu12 = n12,
      progress = 0, Kspec = Kspec, testtau = testtau
    ),
    .progress = progress
  )
}

#' Simulate points according to a bivariate DPP
#'
#' @param rho1 A numeric scalar specifying the intensity of the 1st DPP (default: 100).
#' @param rho2 A numeric scalar specifying the intensity of the 2nd DPP (default: 100).
#' @param tau A numeric scalar specifying the cross-correlation (default: 0.2).
#' @param alpha1 A numeric scalar specifying the ??? of the 1nd DPP (default: 0.03).
#' @param alpha2 A numeric scalar specifying the ??? of the 2nd DPP (default: 0.03).
#' @param alpha12 A numeric scalar specifying the cross-??? (default: 0.05).
#' @param nu1 A numeric scalar specifying the ??? of the first DPP (default: 10).
#' @param nu2 A numeric scalar specifying the ??? of the 2nd DPP (default: 10).
#' @param nu12 A numeric scalar specifying the cross-??? (default: 10).
#' @param progress A numeric specifying whether a progess bar should be displayed (default: 0 for none).
#' @param Kspec A function specifying the kernel to be used (default: \code{Kspecbessel}).
#' @param testtau A function specifying the upper bound for the cross-correlation (default: \code{testtaubessel}).
#'
#' @return A \code{\link[spatstat]{ppp}} object containing the simulated point pattern.
#' @export
#'
#' @examples
#' pp <- rbidpp(progress = 1, Kspec = "Kspecbessel", testtau = "testtaubessel")
rbidpp <- function(
  rho1 = 100, rho2 = 100, tau = 0.2,
  alpha1 = 0.03, alpha2 = 0.03, alpha12 = 0.05,
  nu1 = 10, nu2 = 10, nu12 = 10,
  progress = 0,
  Kspec="Kspecmatern",
  testtau="testtaumatern") {
  testtau<-get(testtau)
  if(!testtau(tau,rho1,rho2,alpha1,alpha2,alpha12)){stop('invalid value for tau')}

  #step1
  expnum=rho1+rho2
  if(Kspec=="Kspecbessel"){prec0<-0.95}
  else{prec0 <- 0.99}
  Kspec<-get(Kspec)

  #prec0 <- 0.95 #for bessel
  trunc <- 1
  prec <- 0
  maxtrunc=1000
  sumdiag=function(r,K,...){sum(diag(K(r,...)))}
  while (prec <= prec0 & (2 * trunc) <= maxtrunc) {
    trunc <- 2 * trunc
    index1a <- c(rep(0, trunc), 1:trunc)
    index2a <- c(1:trunc, rep(0, trunc))
    index1 <- rep(1:trunc, trunc)
    index2 <- rep(1:trunc, each = trunc)
    eigo <- sumdiag(0,Kspec,rho1,rho2,alpha1,alpha2,alpha12,tau)
    eiga <- sapply(sqrt((index1a)^2 + (index2a)^2),sumdiag,Kspec,rho1,rho2,alpha1,alpha2,alpha12,tau)
    eig <- sapply(sqrt((index1)^2 + (index2)^2),sumdiag,Kspec,rho1,rho2,alpha1,alpha2,alpha12,tau)
    prec <- (eigo + 2 * sum(eiga) + 4 * sum(eig))/expnum
  }
  N=max(index1a)


  #step2
  if(progress>0) cat("Spectral decomposition...")
  M=2
  kk=as.matrix(expand.grid(-N:N,-N:N))
  Vfull=NULL
  eigenvalues=NULL
  for(i in 1:nrow(kk)){
    #print(i)
    tmp<-eigen(Kspec(sqrt(kk[i,1]^2+kk[i,2]^2),rho1,rho2,alpha1,alpha2,alpha12,tau))
    eigenvalues=c(eigenvalues,tmp$values)
    Vfull=cbind(Vfull,tmp$vector)
  }
  #dim(Vfull)
  #length(eigenvalues)
  kkfull=kronecker(kk,rep(1,M))
  if(progress>0) cat(" Done.\n")


  #step3
  tmp=rbinom(nrow(kkfull),1,eigenvalues)
  index <-which(tmp==1)
  kkindex=kkfull[index,]
  dim(kkindex)
  V=Vfull[,index]
  #dim(V)
  rm(kkfull)
  rm(Vfull)
  gc()


  #step4 : cf the functions below
  X <- rdpppmulti(kkindex,V,progress=progress)
  return(X)
}




## Generates an empty point pattern
emptyppx <- function(W, simplify = TRUE){
  W <- spatstat::as.boxx(W)
  r <- W$ranges
  d <- ncol(r)
  if(simplify){
    if(d==2)
      return(spatstat::ppp(numeric(0), numeric(0), window = spatstat::as.owin(W)))
    if(d==3)
      return(spatstat::pp3(numeric(0), numeric(0), numeric(0), W))
  }
  rslt <- replicate(d, numeric(0), simplify=FALSE)
  names(rslt) <- paste("x",1:d,sep="")
  rslt <- as.data.frame(rslt)
  return(ppx(rslt, domain = W, coord.type= rep("spatial", d)))
}



##Generates a multidimensional projection DPP  : adaptation from rdppp of spatstat (all changed lines are indicated by -fred)
#rdppp <- function(index, basis = "fourierbasis", window = boxx(rep(list(0:1), ncol(index))), -fred
#                  reject_max = 1e4, progress = 0, debug = FALSE, given = NULL, given_max_volume = 0.5, ...) -fred
rdpppmulti <- function(index, V=NULL, basis = "fourierbasis", window = spatstat::boxx(rep(list(0:1), ncol(index))),
                       reject_max = 1e4, progress = 0, debug = FALSE, given = NULL, given_max_volume = 0.5, ...) #-fred
{
  ##Check if really multi -fred
  if(is.null(V)){return(rdppp(index,basis = "fourierbasis", window = spatstat::boxx(rep(list(0:1),ncol(index))), reject_max = 1e4, progress = 0, debug = FALSE, given = NULL, given_max_volume = 0.5, ...))} #-fred
  if(is.matrix(V) & debug){warning(paste(sQuote("debug"),"is not available for multidimensional DPPs"));debug=FALSE} #-fred
  if(is.matrix(V) & !is.null(given)){warning(paste(sQuote("given"),"is not available for multidimensional DPPs"));given=NULL} #-fred

  ## Check arguments:
  if (!(is.logical(debug)))
    stop(paste(sQuote("debug"), "must be TRUE or FALSE"))
  if (!is.numeric(reject_max)||reject_max<=1)
    stop(paste(sQuote("reject_max"), "must be a numeric greater than 1"))
  if (!is.numeric(progress)||reject_max<1)
    stop(paste(sQuote("progress"), "must be a numeric greater than or equal to 1"))
  index <- as.matrix(index)
  d <- ncol(index)
  window <- spatstat::as.boxx(window)
  ranges <- window$ranges
  boxlengths <- as.numeric(ranges[2L, ] - ranges[1L, ])
  if(ncol(ranges)!=d)
    stop("The dimension differs from the number of columns in index")
  if(basis != "fourierbasis"){
    warning("Non Fourier basis probably doesn't work correctly! Fourier is
            assumed for bounds in rejection sampling.")
    userbasis <- get(basis)
    if (!(is.function(userbasis)))
      stop(paste(sQuote("basis"), "must be a function"))
    tmp <- userbasis(ranges[1,,drop=FALSE], index, window)
    if (!(is.numeric(tmp) || is.complex(tmp)))
      stop(paste("Output of", sQuote("basis"), "must be numeric or complex"))
    basis <- function(x, k, boxlengths){
      userbasis(x, k, spatstat::boxx(lapply(boxlengths, function(x) list(c(0,x)))))
    }
  } else{
    basis <- spatstat::fourierbasisraw
  }

  ## Number of points to simulate:
  n <- nrow(index)
  ## Number of types (=marks) -fred
  M=nrow(V) #-fred

  ## Resolve `given` for pseudo conditional simulation
  if(!is.null(given)){
    # Make sure `given` is a list of point patterns
    if(spatstat::is.ppp(given) || spatstat::is.pp3(given) || spatstat::is.ppx(given)){
      given <- list(given)
    }
    stopifnot(all(sapply(given, function(x){ spatstat::is.ppp(x) || spatstat::is.pp3(x) || spatstat::is.ppx(x)})))
    # Check that the window (or its boundingbox) is inside the simulation window
    Wgiven <- lapply(given, function(x) spatstat::as.boxx(domain(x)))
    stopifnot(all(sapply(Wgiven,
                         function(w){
                           all(w$ranges[1,] >= ranges[1,]) && all(w$ranges[2,] <= ranges[2,])
                         })))
    stopifnot(sum(sapply(Wgiven, volume))<given_max_volume)
    # Resolve number of given points and extract coordinates
    ngiven <- sum(sapply(given, npoints))
    stopifnot(ngiven <= n)
    if(ngiven == n) return(given)
    coordsgiven <- lapply(given, function(x) as.matrix(coords(x)))
    coordsgiven <- Reduce(rbind, coordsgiven)
  }

  ## Return empty point pattern if n=0:
  empty <- emptyppx(window)
  if (n==0)
    return(empty)

  ## Initialize debug info:
  if(debug){
    debugList = replicate(n, list(old=empty, accepted=empty, rejected=empty, index=index), simplify=FALSE)
  }




  # Matrix of coordinates:
  x <- matrix(0,n,d)
  colnames(x) <- paste("x",1:d,sep="")
  x[1,] <- runif(d,as.numeric(ranges[1,]),as.numeric(ranges[2,]))
  if(!is.null(given)){
    x[1,] <- coordsgiven[1,,drop=FALSE]
  }
  type <- rep(0,n) #-fred
  probs<-apply(V^2,1,mean) #-fred
  type[1]<-sample(1:M,1,prob=probs) #-fred

  # Debug info:
  if(debug){
    debugList[[1]]=list(old=empty, accepted=ppx(x[1,,drop=FALSE],window,simplify=TRUE), rejected=empty, index=index, estar=rep(1/n,n))
  }

  if (n==1)
    return(spatstat::ppp(x[1,1],x[1,2], window = spatstat::as.owin(window),marks=as.factor(type[1]))) #-fred
  #return(ppx(x, window, simplify = TRUE))-fred

  # First vector of basis-functions evaluated at first point:
  #v <- basis(x[1,,drop=FALSE],index,boxlengths) -fred
  v <- V[type[1],]*basis(x[1, , drop = FALSE], index,boxlengths) #-fred

  ## Record normalized version in the Gram-Schmidt matrices:
  e <- v/sqrt(sum(abs(v)^2))
  estar <- Conj(e)
  if(progress>0)
    cat(paste("Simulating", n, "points:\n"))

  ## Main for loop over number of points:
  for(i in (n-1):1){
    ## Print progress:
    if(progress>0)
      spatstat::progressreport(n-i, n, every=progress)
    ## Aux. variable to count number of rejection steps:
    tries <- 0
    # Debug info:
    if(debug){
      rejected <- matrix(NA,reject_max,d)
    }
    repeat{
      ## Proposed point:
      newtype <- sample(1:M,1,prob=probs) #-fred
      newx <- matrix(runif(d,as.numeric(ranges[1,]),as.numeric(ranges[2,])),ncol=d)
      if(!is.null(given)){
        if(i>(n-ngiven)){
          newx <- coordsgiven[n-i+1,,drop=FALSE]
        } else{
          while(any(sapply(Wgiven, function(w) spatstat::inside.boxx(as.hyperframe(newx), w = w))))
            newx <- matrix(runif(d,as.numeric(ranges[1,]),as.numeric(ranges[2,])),ncol=d)
        }
      }
      ## Basis functions eval. at proposed point:
      #v <- as.vector(basis(newx, index, boxlengths)) -fred
      v <- as.vector(V[newtype,]*basis(newx, index, boxlengths)) #-fred
      ## Vector of projection weights (has length n-i)
      wei <- t(v)%*%estar
      if(!is.null(given) && i>(n-ngiven)){
        break
      }
      ## Accept probability:
      # tmp <- prod(ranges[2,]-ranges[1,])/n*(sum(abs(v)^2)-sum(abs(wei)^2))
      # tmp <- 1-prod(ranges[2,]-ranges[1,])/n*(sum(abs(wei)^2)) -fred
      tmp <- prod(ranges[2, ] - ranges[1, ])/n * (sum(abs(v)^2) - sum(abs(wei)^2)) / probs[newtype] #-fred

      ## If proposal is accepted the loop is broken:
      if(runif(1)<as.numeric(abs(tmp))){
        break
      }
      ## If rejected, check that we have not tried too many times:
      if(tries>reject_max){
        stop(paste("Rejection sampling failed reject_max =",reject_max,"times in a row"))
      }
      ## Increase the count of rejection steps:
      tries <- tries+1
      # Debug info:
      if(debug){
        rejected[tries,] <- newx
      }
    } ## END OF REJECTION LOOP

    # Record the accepted point:
    x[n-i+1,] <- newx
    type[n - i + 1] <- newtype #-fred


    # Debug info:
    if(debug){
      if(tries==0){
        rej <- empty
      } else{
        rej <- ppx(rejected[1:tries,,drop=FALSE],window, simplify=TRUE)
      }
      debugList[[n-i+1]] = list(
        old=ppx(x[1:(n-i),,drop=FALSE],window, simplify=TRUE),
        accepted=ppx(newx,window,simplify=TRUE),
        rejected=rej, index=index, estar = estar)
    }

    ## If it is the last point exit the main loop:
    if(i==1){break}

    ## Calculate orthogonal vector for Gram-Schmidt procedure:
    # w <- v - rowSums(matrix(wei,n,n-i,byrow=TRUE)*e[,1:(n-i)])
    w <- v - colSums(t(e)*as.vector(wei))
    ## Record normalized version in the Gram-Schmidt matrices:
    enew <- w/sqrt(sum(abs(w)^2))
    e <- cbind(e, enew)
    estar <- cbind(estar,Conj(enew))
  } ## END OF MAIN FOR LOOP
  # Save points as point pattern:
  #X <- ppx(x, window, simplify = TRUE)  -fred
  X <- spatstat::ppp(x[,1],x[,2], window = spatstat::as.owin(window), marks=as.factor(type)) #-fred

  # Debug info:
  if(debug){
    attr(X, "dpp") <- list(debug=debugList)
  }
  if(progress>0)
    cat(" Done!\n")
  return(X)
}
