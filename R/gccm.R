.gccm_sf_method = \(data, cause, effect, libsizes, E = 3, tau = 1, k = E+2, theta = 1, algorithm = "simplex", lib = NULL, pred = NULL,
                    nb = NULL,threads = detectThreads(),parallel.level = "low",bidirectional = TRUE,detrend = TRUE,progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  tau = .check_inputelementnum(tau,2)
  k = .check_inputelementnum(k,2)
  pl = .check_parallellevel(parallel.level)
  .varname = .internal_varname()
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  data = sf::st_drop_geometry(data)
  data = data[,varname]
  names(data) = .varname

  if (detrend){
    data = .internal_detrend(data,.varname,coords)
  }
  cause = data[,"cause",drop = TRUE]
  effect = data[,"effect",drop = TRUE]

  if (is.null(lib)) lib = .internal_library(data)
  if (is.null(pred)) pred = lib

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCCM4Lattice(cause,effect,nb,libsizes,lib,pred,E[1],tau[1],k[1],simplex,theta,threads,pl,progressbar)
  }
  y_xmap_x = RcppGCCM4Lattice(effect,cause,nb,libsizes,lib,pred,E[2],tau[2],k[2],simplex,theta,threads,pl,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

.gccm_spatraster_method = \(data, cause, effect, libsizes, E = 3, tau = 1, k = E+2, theta = 1, algorithm = "simplex", lib = NULL, pred = NULL,
                            threads = detectThreads(), parallel.level = "low", bidirectional = TRUE, detrend = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  tau = .check_inputelementnum(tau,2)
  k = .check_inputelementnum(k,2)
  pl = .check_parallellevel(parallel.level)
  libsizes = as.matrix(libsizes)
  .varname = .internal_varname()
  data = data[[varname]]
  names(data) = .varname

  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (detrend){
    dtf = .internal_detrend(dtf,.varname)
  }
  causemat = matrix(dtf[,"cause"],nrow = terra::nrow(data),byrow = TRUE)
  effectmat = matrix(dtf[,"effect"],nrow = terra::nrow(data),byrow = TRUE)

  if (is.null(lib)) lib = .internal_library(dtf,TRUE)
  if (is.null(pred)) pred = lib

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCCM4Grid(causemat,effectmat,libsizes,lib,pred,E[1],tau[1],k[1],simplex,theta,threads,pl,progressbar)
  }
  y_xmap_x = RcppGCCM4Grid(effectmat,causemat,libsizes,lib,pred,E[2],tau[2],k[2],simplex,theta,threads,pl,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' geographical convergent cross mapping
#'
#' @param data observation data.
#' @param cause name of causal variable.
#' @param effect name of effect variable.
#' @param libsizes number of spatial units used in prediction.
#' @param E (optional) embedding dimensions.
#' @param tau (optional) step of spatial lags.
#' @param k (optional) number of nearest neighbors used in prediction.
#' @param theta (optional) weighting parameter for distances, useful when `algorithm` is `smap`.
#' @param algorithm (optional) prediction algorithm.
#' @param lib (optional) libraries indices.
#' @param pred (optional) predictions indices.
#' @param nb (optional) neighbours list.
#' @param threads (optional) number of threads to use.
#' @param parallel.level (optional) level of parallelism, `low` or `high`.
#' @param bidirectional (optional) whether to examine bidirectional causality.
#' @param detrend (optional) whether to remove the linear trend.
#' @param progressbar (optional) whether to show the progress bar.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping results}
#' \item{\code{varname}}{names of causal and effect variable}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name gccm
#' @references
#' Gao, B., Yang, J., Chen, Z. et al. Causal inference from cross-sectional earth system data with geographical convergent cross mapping. Nat Commun 14, 5875 (2023).
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' g = gccm(columbus,"hoval","crime",libsizes = seq(5,45,5),E = 6)
#' g
#' plot(g, ylimits = c(0,0.85))
#' }
methods::setMethod("gccm", "sf", .gccm_sf_method)

#' @rdname gccm
methods::setMethod("gccm", "SpatRaster", .gccm_spatraster_method)
