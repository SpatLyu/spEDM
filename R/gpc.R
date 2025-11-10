.gpc_sf_method = \(data, cause, effect, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, dist.metric = "L2", zero.tolerance = E+2,
                   relative = TRUE, weighted = TRUE, na.rm = TRUE, threads = detectThreads(), detrend = FALSE, bidirectional = TRUE, nb = NULL){
  varname = .check_character(cause, effect)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  cause = .uni_lattice(data,cause,detrend)
  effect = .uni_lattice(data,effect,detrend)
  dist.metric = .check_distmetric(dist.metric)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)))
  if (is.null(pred)) pred = lib
  return(.run_gpc(cause, effect, E, k, tau, style, lib, pred,
                  dist.metric, zero.tolerance,
                  relative, weighted, na.rm, threads,
                  bidirectional, varname, nb))
}

.gpc_spatraster_method = \(data, cause, effect, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, dist.metric = "L2", zero.tolerance = E+2,
                           relative = TRUE, weighted = TRUE, na.rm = TRUE, threads = detectThreads(), detrend = FALSE, bidirectional = TRUE, grid.coord = TRUE){
  varname = .check_character(cause, effect)
  cause = .uni_grid(data,cause,detrend,grid.coord)
  effect = .uni_grid(data,effect,detrend,grid.coord)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  return(.run_gpc(cause, effect, E, k, tau, style, lib, pred,
                  dist.metric, zero.tolerance,
                  relative, weighted, na.rm, threads,
                  bidirectional, varname))
}

#' geographical pattern causality
#'
#' @param data observation data.
#' @param cause name of causal variable.
#' @param effect name of effect variable.
#' @param E (optional) embedding dimension.
#' @param k (optional) number of nearest neighbors used in cross mapping.
#' @param tau (optional) step of spatial lag.
#' @param style (optional) embedding style (`0` includes current state, `1` excludes it).
#' @param lib (optional) libraries indices (input needed: `vector` - spatial vector, `matrix` - spatial raster).
#' @param pred (optional) predictions indices (input requirement same as `lib`).
#' @param dist.metric (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).
#' @param zero.tolerance (optional) maximum number of zeros tolerated in signature space.
#' @param relative (optional) whether to calculate relative changes in embeddings.
#' @param weighted (optional) whether to weight causal strength.
#' @param na.rm (optional) whether to remove `NA` samples in symbolic pattern generation.
#' @param threads (optional) number of threads to use.
#' @param detrend (optional) whether to remove the linear trend.
#' @param bidirectional (optional) whether to examine bidirectional causality.
#' @param nb (optional) neighbours list.
#'
#' @return A list
#' \describe{
#' \item{\code{causality}}{per-sample causality statistics}
#' \item{\code{summary}}{overall causal strength}
#' \item{\code{pattern}}{pairwise pattern relationships}
#' \item{\code{varname}}{names of causal and effect variable}
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
#' }
methods::setMethod("gpc", "sf", .gpc_sf_method)

#' @rdname gpc
#' @param grid.coord (optional) whether to detrend using cell center coordinates (`TRUE`) or row/column numbers (`FALSE`).
methods::setMethod("gpc", "SpatRaster", .gpc_spatraster_method)
