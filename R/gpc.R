.gpc_sf_method = \(data, cause, effect, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99, random = TRUE, seed = 42, dist.metric = "L2", zero.tolerance = k,
                   relative = TRUE, weighted = TRUE, na.rm = TRUE, threads = detectThreads(), detrend = FALSE, parallel.level = "low", bidirectional = TRUE, progressbar = TRUE, nb = NULL){
  varname = .check_character(cause, effect)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  cause = .uni_lattice(data,cause,detrend)
  effect = .uni_lattice(data,effect,detrend)
  dist.metric = .check_distmetric(dist.metric)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)))
  if (is.null(pred)) pred = lib
  return(.run_gpc(cause, effect, E, k, tau, style, lib, pred, dist.metric,
                  zero.tolerance, relative, weighted, na.rm, threads, bidirectional,
                  varname, nb, libsizes, boot, random, seed, parallel.level, progressbar))
}

.gpc_spatraster_method = \(data, cause, effect, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99, random = TRUE, seed = 42, dist.metric = "L2", zero.tolerance = k,
                           relative = TRUE, weighted = TRUE, na.rm = TRUE, threads = detectThreads(), detrend = FALSE, parallel.level = "low", bidirectional = TRUE, progressbar = TRUE, grid.coord = TRUE){
  varname = .check_character(cause, effect)
  cause = .uni_grid(data,cause,detrend,grid.coord)
  effect = .uni_grid(data,effect,detrend,grid.coord)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  return(.run_gpc(cause, effect, E, k, tau, style, lib, pred, dist.metric,
                  zero.tolerance, relative, weighted, na.rm, threads, bidirectional,
                  varname, NULL, libsizes, boot, random, seed, parallel.level, progressbar))
}

#' geographical pattern causality
#'
#' @inheritParams gcmc
#' @param boot (optional) number of bootstraps to perform.
#' @param seed (optional) random seed.
#' @param random (optional) whether to use random sampling.
#' @param zero.tolerance (optional) maximum number of zeros tolerated in signature space.
#' @param relative (optional) whether to calculate relative changes in embeddings.
#' @param weighted (optional) whether to weight causal strength.
#' @param na.rm (optional) whether to remove `NA` samples in symbolic pattern generation.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping results (only present if `libsizes` is not `NULL`)}
#' \item{\code{causality}}{per-sample causality statistics(present if `libsizes` is `NULL`)}
#' \item{\code{summary}}{overall causal strength(present if `libsizes` is `NULL`)}
#' \item{\code{pattern}}{pairwise pattern relationships(present if `libsizes` is `NULL`)}
#' \item{\code{varname}}{names of causal and effect variable}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name gpc
#' @aliases gpc,sf-method
#' @references
#' Zhang, Z., Wang, J., 2025. A model to identify causality for geographic patterns. International Journal of Geographical Information Science 1–21.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' gpc(columbus,"hoval","crime", E = 3)
#'
#' # convergence diagnostics
#' g = gpc(columbus,"hoval","crime",libsizes = seq(5,45,5),E = 3)
#' plot(g)
#' }
methods::setMethod("gpc", "sf", .gpc_sf_method)

#' @rdname gpc
#' @param grid.coord (optional) whether to detrend using cell center coordinates (`TRUE`) or row/column numbers (`FALSE`).
methods::setMethod("gpc", "SpatRaster", .gpc_spatraster_method)
