.sc_sf_method = \(data, cause, effect, k, block = 3, boot = 399, seed = 42, base = 2, lib = NULL, pred = NULL,
                  nb = NULL, threads = detectThreads(), detrend = TRUE, normalize = FALSE, progressbar = FALSE){
  varname = .check_character(cause, effect)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  block = RcppDivideLattice(nb,block)
  cause = .uni_lattice(data,cause,detrend)
  effect = .uni_lattice(data,effect,detrend)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)))
  if (is.null(pred)) pred = lib
  return(.bind_sc(RcppSGC4Lattice(cause,effect,nb,lib,pred,block,k,threads,boot,base,seed,TRUE,normalize,progressbar),varname))
}

.sc_spatraster_method = \(data, cause, effect, k, block = 3, boot = 399, seed = 42, base = 2, lib = NULL, pred = NULL,
                          threads = detectThreads(), detrend = TRUE, normalize = FALSE, progressbar = FALSE){
  varname = .check_character(cause, effect)
  cause = .uni_grid(data,cause,detrend)
  effect = .uni_grid(data,effect,detrend)
  block = matrix(RcppDivideGrid(effect,block),ncol = 1)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  return(.bind_sc(RcppSGC4Grid(cause,effect,lib,pred,block,k,threads,boot,base,seed,TRUE,normalize,progressbar),varname))
}

#' spatial causality test
#'
#' @param data observation data.
#' @param cause name of causal variable.
#' @param effect name of effect variable.
#' @param k (optional) number of nearest neighbors used in symbolization.
#' @param block (optional) number of blocks used in spatial block bootstrap.
#' @param boot (optional) number of bootstraps to perform.
#' @param seed (optional) random seed.
#' @param base (optional) logarithm base.
#' @param lib (optional) libraries indices.
#' @param pred (optional) predictions indices.
#' @param nb (optional) neighbours list.
#' @param threads (optional) number of threads to use.
#' @param detrend (optional) whether to remove the linear trend.
#' @param normalize (optional) whether to normalize the result.
#' @param progressbar (optional) whether to show the progress bar.
#'
#' @return A list
#' \describe{
#' \item{\code{sc}}{statistic for spatial causality}
#' \item{\code{varname}}{names of causal and effect variable}
#' }
#' @export
#' @name sc.test
#' @aliases sc.test,sf-method
#' @references
#' Herrera, M., Mur, J., & Ruiz, M. (2016). Detecting causal relationships between spatial processes. Papers in Regional Science, 95(3), 577–595.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' sc.test(columbus,"hoval","crime", k = 15)
#' }
methods::setMethod("sc.test", "sf", .sc_sf_method)

#' @rdname sc.test
methods::setMethod("sc.test", "SpatRaster", .sc_spatraster_method)
