methods::setGeneric("sc.test", function(data, ...) standardGeneric("sc.test"))

.sct_sf_method = \(data, cause, effect, k, block = 3, boot = 399, seed = 42, base = 2, lib = NULL, pred = NULL,
                   nb = NULL, threads = detectThreads(), symbolize = TRUE, normalize = FALSE, progressbar = FALSE){
  varname = .check_character(cause, effect)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  block = RcppDivideLattice(nb,block)
  cause = .uni_lattice(data,cause,FALSE)
  effect = .uni_lattice(data,effect,FALSE)
  return(.bind_sct(RcppSCT4Lattice(cause,effect,nb,lib,pred,block,k,threads,boot,base,seed,symbolize,normalize,progressbar),varname))
}

.sct_spatraster_method = \(data, cause, effect, k, block = 3, boot = 399, seed = 42, base = 2, lib = NULL, pred = NULL,
                           threads = detectThreads(),symbolize = TRUE,normalize = FALSE,progressbar = FALSE){
  varname = .check_character(cause, effect)
  cause = .uni_grid(data,cause,FALSE)
  effect = .uni_grid(data,effect,FALSE)
  block = matrix(RcppDivideGrid(effect,block),ncol = 1)
  if (is.null(lib)) lib = which(!is.na(mat), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  return(.bind_sct(RcppSCT4Grid(cause,effect,lib,pred,block,k,threads,boot,base,seed,symbolize,normalize,progressbar),varname))
}

#' spatial (granger) causality test
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
#' @param symbolize (optional) Whether to apply the symbolic map process.
#' @param normalize (optional) Whether to normalize the result to `[-1, 1]`.
#' @param progressbar (optional) Whether to show the progress bar.
#'
#' @return A list
#' \describe{
#' \item{\code{sc}}{statistic for spatial causality}
#' \item{\code{varname}}{names of causal and effect variable}
#' }
#' @export
#' @name sc.test
#' @rdname sc.test
#' @aliases sc.test,sf-method
#' @references
#' Herrera, M., Mur, J., & Ruiz, M. (2016). Detecting causal relationships between spatial processes. Papers in Regional Science, 95(3), 577–595.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' sc.test(columbus,"hoval","crime", k = 15)
#' }
methods::setMethod("sc.test", "sf", .sct_sf_method)

#' @rdname sc.test
methods::setMethod("sc.test", "SpatRaster", .sct_spatraster_method)
