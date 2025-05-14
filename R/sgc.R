methods::setGeneric("sgc", function(data, ...) standardGeneric("sgc"))

.sgc_sf_method = \(data, cause, effect, k, block = 3, boot = 399, seed = 42, base = 2, lib = NULL, pred = NULL,
                   nb = NULL, threads = detectThreads(), trend.rm = TRUE, normalize = FALSE, progressbar = FALSE){
  varname = .check_character(cause, effect)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  block = RcppDivideLattice(nb,block)
  cause = .uni_lattice(data,cause,trend.rm)
  effect = .uni_lattice(data,effect,trend.rm)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)))
  if (is.null(pred)) pred = lib
  return(.bind_sc(RcppSGC4Lattice(cause,effect,nb,lib,pred,block,k,threads,boot,base,seed,TRUE,normalize,progressbar),varname))
}

.sgc_spatraster_method = \(data, cause, effect, k, block = 3, boot = 399, seed = 42, base = 2, lib = NULL, pred = NULL,
                           threads = detectThreads(), trend.rm = TRUE, normalize = FALSE, progressbar = FALSE){
  varname = .check_character(cause, effect)
  cause = .uni_grid(data,cause,trend.rm)
  effect = .uni_grid(data,effect,trend.rm)
  block = matrix(RcppDivideGrid(effect,block),ncol = 1)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  return(.bind_sc(RcppSGC4Grid(cause,effect,lib,pred,block,k,threads,boot,base,seed,TRUE,normalize,progressbar),varname))
}

#' spatial granger causality test
#'
#' @param data The observation data.
#' @param cause Name of causal variable.
#' @param effect Name of effect variable.
#' @param k (optional) Number of nearest neighbors used in symbolization.
#' @param block (optional) Number of blocks used in spatial block bootstrap.
#' @param boot (optional) Number of bootstraps to perform.
#' @param seed (optional) The random seed.
#' @param base (optional) Base of the logarithm.
#' @param lib (optional) Libraries indices.
#' @param pred (optional) Predictions indices.
#' @param nb (optional) The neighbours list.
#' @param threads (optional) Number of threads.
#' @param trend.rm (optional) .
#' @param normalize (optional) Whether to normalize the result.
#' @param progressbar (optional) Whether to show the progress bar.
#'
#' @return A list
#' \describe{
#' \item{\code{sc}}{statistic for spatial causality}
#' \item{\code{varname}}{names of causal and effect variable}
#' }
#' @export
#' @name sgc
#' @rdname sgc
#' @aliases sgc,sf-method
#' @references
#' Herrera, M., Mur, J., & Ruiz, M. (2016). Detecting causal relationships between spatial processes. Papers in Regional Science, 95(3), 577–595.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' sgc(columbus,"hoval","crime", k = 15)
#' }
methods::setMethod("sgc", "sf", .sgc_sf_method)

#' @rdname sgc
methods::setMethod("sgc", "SpatRaster", .sgc_spatraster_method)
