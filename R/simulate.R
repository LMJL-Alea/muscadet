#' Random generation of point patterns from a cross-type DPP
#'
#' Simulates one or more planar point patterns following a given multi-mark DPP.
#'
#' @param n An integer value specifying the sample size. Defaults to `1L`.
#' @param rho1 A numeric value specifying the intensity of the 1st marginal DPP.
#'   Defaults to `100L`.
#' @param rho2 A numeric value specifying the intensity of the 2nd marginal DPP.
#'   Defaults to `100L`.
#' @param tau A numeric value specifying the correlation between marks. Defaults
#'   to `0.2`.
#' @param alpha1 A numeric value specifying the repulsion rate of the 1st DPP.
#'   Defaults to `0.03`.
#' @param alpha2 A numeric value specifying the repulsion rate of the 2nd DPP.
#'   Defaults to `0.03`.
#' @param alpha12 A numeric value specifying the repulsion rate between marks.
#'   Defaults to `0.05`.
#' @param L An integer value specifying the window size of the point pattern.
#'   Defaults to `1L`.
#' @param d An integer value specifying the dimension of the DPP. Defaults to
#'   `2L`.
#' @param model A string specifying the model to be used. Choices are `"Gauss"`
#'   or `"Bessel"`. Defaults to `"Gauss"`.
#' @param max_trunc An integer value specifying the maximum number of basis
#'   functions in the Fourier decomposition. Defaults to `1000L`.
#'
#' @return An object of class [spatstat.geom::ppp] storing the simulated point
#'   pattern.
#'
#' @keywords internal
rbidpp <- function(n = 1L,
                   rho1 = 100, rho2 = 100,
                   alpha1 = 0.03, alpha2 = 0.03,
                   alpha12 = 0.035, tau = 0.5,
                   L = 1, d = 2L, model = c("Gauss", "Bessel"),
                   max_trunc = 1000L)
{
  pb <- progressr::progressor(along = 1:n)
  model <- rlang::arg_match(model)
  furrr::future_map(1:n, ~ {
    pb()
    .rbidpp_single(
      rho1 = rho1, rho2 = rho2,
      alpha1 = alpha1, alpha2 = alpha2,
      alpha12 = alpha12, tau = tau,
      L = L, d = d, model = model,
      maxtrunc = max_trunc
    )
  }, .options = furrr::furrr_options(seed = TRUE))
}

sumdiag <- function(r, K, L, ...) {
  sum(diag(K(r / L, ...)))
}

.rbidpp_single <- function(rho1 = 100, rho2 = 100,
                           alpha1 = 0.03, alpha2 = 0.03,
                           alpha12 = 0.035, tau = 0.5,
                           L = 1, d = 2, model = c("Gauss", "Bessel"),
                           maxtrunc = 1000) {
  if (!check_parameter_set(rho1, rho2, alpha1, alpha2, alpha12, tau, d, model))
    cli::cli_abort("The set of parameters is not valid for a 2-mark DPP.")

  model <- rlang::arg_match(model)

  # Step1
  expnum <- (rho1 + rho2) * L ^ 2
  prec0 <- switch(
    model,
    Gauss = 0.99,
    Bessel = 0.95
  )

  trunc <- 1
  prec <- 0
  while (prec <= prec0 & (2 * trunc) <= maxtrunc) {
    trunc <- 2 * trunc
    index1a <- c(rep(0, trunc), 1:trunc)
    index2a <- c(1:trunc, rep(0, trunc))
    index1 <- rep(1:trunc, trunc)
    index2 <- rep(1:trunc, each = trunc)
    eigo <- sumdiag(
      r = 0, K = get_khat_matrix, L = L,
      rho1 = rho1, rho2 = rho2,
      alpha1 = alpha1, alpha2 = alpha2,
      alpha12 = alpha12, tau = tau,
      d = d, model = model
    )
    eiga <- purrr::map_dbl(
      .x = sqrt((index1a) ^ 2 + (index2a) ^ 2),
      .f = sumdiag,
      K = get_khat_matrix, L = L,
      rho1 = rho1, rho2 = rho2,
      alpha1 = alpha1, alpha2 = alpha2,
      alpha12 = alpha12, tau = tau,
      d = d, model = model
    )
    eig <- purrr::map_dbl(
      .x = sqrt((index1) ^ 2 + (index2) ^ 2),
      .f = sumdiag,
      K = get_khat_matrix, L = L,
      rho1 = rho1, rho2 = rho2,
      alpha1 = alpha1, alpha2 = alpha2,
      alpha12 = alpha12, tau = tau,
      d = d, model = model
    )
    prec <- (eigo + 2 * sum(eiga) + 4 * sum(eig)) / expnum
  }

  N <- max(index1a)

  # Step2
  l <- rbidpp_impl(
    N = N, L = L,
    rho1 = rho1, rho2 = rho2,
    alpha1 = alpha1, alpha2 = alpha2,
    alpha12 = alpha12, tau = tau,
    model = model, nbThreads = 1
  )
  kkindex <- l$kkindex
  V <- l$V

  ordering_idx <- kkindex |>
    `colnames<-`(c("x", "y")) |>
    tibble::as_tibble() |>
    dplyr::mutate(norm2 = .data$x^2 + .data$y^2, index = 1:dplyr::n()) |>
    dplyr::arrange(.data$norm2, .data$x, .data$y) |>
    dplyr::pull(.data$index)
  kkindex <- kkindex[ordering_idx, ]
  V <- V[, ordering_idx]

  # Step3 : cf the functions below
  rdpppmulti(
    index = kkindex,
    V = V,
    window = spatstat.geom::boxx(purrr::map(1:ncol(kkindex), ~ c(0, L)))
  )
}

## Generates an empty point pattern
emptyppx <- function(W, simplify = TRUE){
  W <- spatstat.geom::as.boxx(W)
  r <- W$ranges
  d <- ncol(r)
  if(simplify){
    if(d==2)
      return(spatstat.geom::ppp(numeric(0), numeric(0), window = spatstat.geom::as.owin(W)))
    if(d==3)
      return(spatstat.geom::pp3(numeric(0), numeric(0), numeric(0), W))
  }
  rslt <- replicate(d, numeric(0), simplify=FALSE)
  names(rslt) <- paste("x",1:d,sep="")
  rslt <- as.data.frame(rslt)
  return(spatstat.geom::ppx(rslt, domain = W, coord.type= rep("spatial", d)))
}

##Generates a multidimensional projection DPP  : adaptation from rdppp of spatstat (all changed lines are indicated by -fred)
#rdppp <- function(index, basis = "fourierbasis", window = boxx(rep(list(0:1), ncol(index))), -fred
#                  reject_max = 1e4, progress = 0, debug = FALSE, given = NULL, given_max_volume = 0.5, ...) -fred
rdpppmulti <- function(index,
                       V = NULL,
                       basis = "fourierbasis",
                       window = spatstat.geom::boxx(rep(list(0:1), ncol(index))),
                       reject_max = 1e5,
                       progress = 0,
                       debug = FALSE,
                       given = NULL,
                       given_max_volume = 0.5,
                       ...) #-fred
{
  ##Check if really multi -fred
  if (is.null(V)) {
    return(spatstat.model::rdpp(
      index = index,
      basis = "fourierbasis",
      window = spatstat.geom::boxx(rep(list(0:1),ncol(index))),
      reject_max = 1e4,
      progress = 0,
      debug = FALSE,
      given = NULL,
      given_max_volume = 0.5,
      ...))
  } #-fred

  if (is.matrix(V) & debug) {
    warning(paste(sQuote("debug"),"is not available for multidimensional DPPs"))
    debug <- FALSE
  } #-fred

  if (is.matrix(V) & !is.null(given)) {
    warning(paste(sQuote("given"),"is not available for multidimensional DPPs"))
    given <- NULL
  } #-fred

  ## Check arguments:
  if (!(is.logical(debug)))
    stop(paste(sQuote("debug"), "must be TRUE or FALSE"))
  if (!is.numeric(reject_max) || reject_max <= 1)
    stop(paste(sQuote("reject_max"), "must be a numeric greater than 1"))
  if (!is.numeric(progress) || reject_max < 1)
    stop(paste(sQuote("progress"), "must be a numeric greater than or equal to 1"))

  index <- as.matrix(index)
  d <- ncol(index)
  window <- spatstat.geom::as.boxx(window)
  ranges <- window$ranges
  boxlengths <- as.numeric(ranges[2L, ] - ranges[1L, ])

  if (ncol(ranges) != d)
    stop("The dimension differs from the number of columns in index")
  if (basis != "fourierbasis") {
    warning("Non Fourier basis probably doesn't work correctly! Fourier is
            assumed for bounds in rejection sampling.")
    userbasis <- get(basis)
    if (!(is.function(userbasis)))
      stop(paste(sQuote("basis"), "must be a function"))
    tmp <- userbasis(ranges[1, , drop = FALSE], index, window)
    if (!(is.numeric(tmp) || is.complex(tmp)))
      stop(paste("Output of", sQuote("basis"), "must be numeric or complex"))
    basis <- function(x, k, boxlengths) {
      userbasis(
        x,
        k,
        spatstat.geom::boxx(lapply(boxlengths, function(x) list(c(0,x))))
      )
    }
  } else {
    basis <- spatstat.geom::fourierbasisraw
  }

  ## Number of points to simulate:
  n <- nrow(index)
  ## Number of types (=marks) -fred
  M <- nrow(V) #-fred

  ## Resolve `given` for pseudo conditional simulation
  if (!is.null(given)) {
    # Make sure `given` is a list of point patterns
    if (spatstat.geom::is.ppp(given) ||
        spatstat.geom::is.pp3(given) ||
        spatstat.geom::is.ppx(given)) {
      given <- list(given)
    }

    stopifnot(all(sapply(given, function(x) {
      spatstat.geom::is.ppp(x) ||
        spatstat.geom::is.pp3(x) ||
        spatstat.geom::is.ppx(x)
    })))

    # Check that the window (or its boundingbox) is inside the simulation window
    Wgiven <- lapply(given, function(x) spatstat.geom::as.boxx(spatstat.geom::domain(x)))
    stopifnot(all(sapply(Wgiven, function(w) {
      all(w$ranges[1, ] >= ranges[1, ]) &&
        all(w$ranges[2, ] <= ranges[2, ])
    })))
    stopifnot(sum(sapply(Wgiven, spatstat.geom::volume)) < given_max_volume)
    # Resolve number of given points and extract coordinates
    ngiven <- sum(sapply(given, spatstat.geom::npoints))
    stopifnot(ngiven <= n)
    if (ngiven == n) return(given)
    coordsgiven <- lapply(given, function(x) as.matrix(spatstat.geom::coords(x)))
    coordsgiven <- Reduce(rbind, coordsgiven)
  }

  ## Return empty point pattern if n=0:
  empty <- emptyppx(window)
  if (n == 0)
    return(empty)

  ## Initialize debug info:
  if (debug) {
    debugList = replicate(
      n = n,
      expr = list(
        old = empty,
        accepted = empty,
        rejected = empty,
        index = index
      ),
      simplify = FALSE
    )
  }

  # Matrix of coordinates:
  x <- matrix(0, n, d)
  colnames(x) <- paste("x", 1:d, sep = "")
  x[1, ] <- stats::runif(d, as.numeric(ranges[1, ]), as.numeric(ranges[2, ]))
  if (!is.null(given))
    x[1, ] <- coordsgiven[1, , drop = FALSE]
  type <- rep(0, n) #-fred
  probs <- apply(V ^ 2, 1, mean) #-fred
  type[1] <- sample(1:M, 1, prob = probs) #-fred

  # Debug info:
  if (debug) {
    debugList[[1]] <- list(
      old = empty,
      accepted = spatstat.geom::ppx(x[1, , drop = FALSE], window, simplify = TRUE),
      rejected = empty,
      index = index,
      estar = rep(1/n, n)
    )
  }

  if (n == 1)
    return(spatstat.geom::ppp(
      x = x[1, 1],
      y = x[1, 2],
      window = spatstat.geom::as.owin(window),
      marks = as.factor(type[1])
    )) #-fred
  #return(ppx(x, window, simplify = TRUE))-fred

  # First vector of basis-functions evaluated at first point:
  #v <- basis(x[1,,drop=FALSE],index,boxlengths) -fred
  v <- V[type[1], ] * basis(x[1, , drop = FALSE], index, boxlengths) #-fred

  ## Record normalized version in the Gram-Schmidt matrices:
  e <- v / sqrt(sum(abs(v) ^ 2))
  estar <- Conj(e)
  if (progress > 0)
    cli::cli_alert_info("Simulating {n} point{?s}")

  ## Main for loop over number of points:

  for (i in (n - 1):1) {
    ## Print progress:
    if (progress > 0)
      spatstat.geom::progressreport(n - i, n, every = progress)
    ## Aux. variable to count number of rejection steps:
    tries <- 0
    # Debug info:
    if(debug){
      rejected <- matrix(NA, reject_max, d)
    }
    repeat {
      ## Proposed point:
      newtype <- sample(1:M, 1, prob = probs) #-fred
      newx <- matrix(
        stats::runif(d, as.numeric(ranges[1, ]), as.numeric(ranges[2, ])),
        ncol = d
      )

      if (!is.null(given)) {
        if (i > (n - ngiven)) {
          newx <- coordsgiven[n - i + 1, , drop = FALSE]
        } else {
          while(any(sapply(Wgiven, function(w) {
            spatstat.geom::inside.boxx(spatstat.geom::as.hyperframe(newx), w = w)
          })))
            newx <- matrix(
              stats::runif(d, as.numeric(ranges[1, ]), as.numeric(ranges[2, ])),
              ncol = d
            )
        }
      }
      ## Basis functions eval. at proposed point:
      #v <- as.vector(basis(newx, index, boxlengths)) -fred
      v <- as.vector(V[newtype,]*basis(newx, index, boxlengths)) #-fred
      ## Vector of projection weights (has length n-i)
      wei <- t(v) %*% estar
      if (!is.null(given) && i > (n - ngiven))
        break

      ## Accept probability:
      # tmp <- prod(ranges[2,]-ranges[1,])/n*(sum(abs(v)^2)-sum(abs(wei)^2))
      # tmp <- 1-prod(ranges[2,]-ranges[1,])/n*(sum(abs(wei)^2)) -fred
      tmp <- prod(ranges[2, ] - ranges[1, ])/n * (sum(abs(v)^2) - sum(abs(wei)^2)) / probs[newtype] #-fred

      ## If proposal is accepted the loop is broken:
      if (stats::runif(1) < as.numeric(abs(tmp)))
        break

      ## If rejected, check that we have not tried too many times:
      if (tries > reject_max)
        cli::cli_abort("Rejection sampling failed {reject_max} times in a row.")
      ## Increase the count of rejection steps:
      tries <- tries + 1
      # Debug info:
      if (debug)
        rejected[tries,] <- newx
    } ## END OF REJECTION LOOP

    # Record the accepted point:
    x[n - i + 1, ] <- newx
    type[n - i + 1] <- newtype #-fred


    # Debug info:
    if (debug) {
      if (tries == 0)
        rej <- empty
      else
        rej <- spatstat.geom::ppx(rejected[1:tries, , drop = FALSE], window, simplify = TRUE)
      debugList[[n - i + 1]] <- list(
        old = spatstat.geom::ppx(x[1:(n - i), ,drop = FALSE], window, simplify = TRUE),
        accepted = spatstat.geom::ppx(newx, window, simplify = TRUE),
        rejected = rej,
        index = index,
        estar = estar
      )
    }

    ## If it is the last point exit the main loop:
    if (i == 1)
      break

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
  X <- spatstat.geom::ppp(
    x = x[, 1],
    y = x[, 2],
    window = spatstat.geom::as.owin(window),
    marks = as.factor(type)
  ) #-fred

  # Debug info:
  if(debug){
    attr(X, "dpp") <- list(debug=debugList)
  }
  if(progress>0)
    cat(" Done!\n")
  return(X)
}
