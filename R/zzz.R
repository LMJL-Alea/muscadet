Rcpp::loadModule("DPP", TRUE)

.onUnload <- function (libpath) {
  library.dynam.unload("mediator", libpath)
}
