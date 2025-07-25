.gcmc_sf_method = \(data, cause, effect, libsizes = NULL, E = 3, tau = 1, k = pmin(E^2), lib = NULL, pred = NULL, nb = NULL,
                    threads = detectThreads(), parallel.level = "low", bidirectional = TRUE, detrend = FALSE, progressbar = TRUE){
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
  if (is.null(libsizes)) libsizes = length(lib)

  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCMC4Lattice(cause,effect,nb,libsizes,lib,pred,E,tau,k[1],0,threads,pl,progressbar)
  }
  y_xmap_x = RcppGCMC4Lattice(effect,cause,nb,libsizes,lib,pred,rev(E),rev(tau),k[2],0,threads,pl,progressbar)

  return(.bind_intersectdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

.gcmc_spatraster_method = \(data, cause, effect, libsizes = NULL, E = 3, tau = 1, k = pmin(E^2), lib = NULL, pred = NULL,
                            threads = detectThreads(), parallel.level = "low", bidirectional = TRUE, detrend = FALSE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  tau = .check_inputelementnum(tau,2)
  k = .check_inputelementnum(k,2)
  pl = .check_parallellevel(parallel.level)
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
  if (is.null(libsizes)) libsizes = matrix(nrow(lib))

  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCMC4Grid(causemat,effectmat,libsizes,lib,pred,E,tau,k[1],0,threads,pl,progressbar)
  }
  y_xmap_x = RcppGCMC4Grid(effectmat,causemat,libsizes,lib,pred,rev(E),rev(tau),k[2],0,threads,pl,progressbar)

  return(.bind_intersectdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' geographical cross mapping cardinality
#'
#' @param data observation data.
#' @param cause name of causal variable.
#' @param effect name of effect variable.
#' @param libsizes (optional) number of spatial units used.
#' @param E (optional) embedding dimensions.
#' @param tau (optional) step of spatial lags.
#' @param k (optional) number of nearest neighbors.
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
#' \item{\code{cs}}{causal strength}
#' \item{\code{varname}}{names of causal and effect variable}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name gcmc
#' @aliases gcmc,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' g = gcmc(columbus,"hoval","crime",E = 7,k = 18)
#' g
#' }
methods::setMethod("gcmc", "sf", .gcmc_sf_method)

#' @rdname gcmc
methods::setMethod("gcmc", "SpatRaster", .gcmc_spatraster_method)
