methods::setGeneric("gccm", function(data, ...) standardGeneric("gccm"))

.gccm_sf_method = \(data, cause, effect, libsizes, E = c(3,3), tau = 1, k = 4, theta = 1, algorithm = "simplex", pred = NULL,
                    nb = NULL, threads = detectThreads(), bidirectional = TRUE, trend.rm = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  k = .check_inputelementnum(k,2)
  tau = .check_inputelementnum(tau,2)
  .varname = .internal_varname()
  lib = 1:nrow(data)
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  if (nrow(data) != length(nb)) stop("Incompatible Data Dimensions!")
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  data = sf::st_drop_geometry(data)
  data = data[,varname]
  names(data) = .varname

  if (trend.rm){
    data = .internal_trend_rm(data,.varname,coords)
  }

  cause = data[,"cause",drop = TRUE]
  effect = data[,"effect",drop = TRUE]

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCCM4Lattice(cause,effect,nb,libsizes,lib,pred,E[1],tau[1],k[1],simplex,theta,threads,progressbar)
  }
  y_xmap_x = RcppGCCM4Lattice(effect,cause,nb,libsizes,lib,pred,E[2],tau[2],k[2],simplex,theta,threads,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

.gccm_spatraster_method = \(data, cause, effect, libsizes, E = c(3,3), tau = 1, k = 4, theta = 1, algorithm = "simplex",
                            pred = NULL, threads = detectThreads(), bidirectional = TRUE, trend.rm = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  k = .check_inputelementnum(k,2)
  tau = .check_inputelementnum(tau,2)
  libsizes = .check_libsizes_gird(libsizes)
  .varname = .internal_varname()
  data = data[[varname]]
  names(data) = .varname

  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (trend.rm){
    dtf = .internal_trend_rm(dtf,.varname)
  }
  causemat = matrix(dtf[,"cause"],nrow = terra::nrow(data),byrow = TRUE)
  effectmat = matrix(dtf[,"effect"],nrow = terra::nrow(data),byrow = TRUE)

  if (is.null(pred)) pred = .internal_predmat(causemat)

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCCM4Grid(causemat,effectmat,libsizes,pred,E[1],tau[1],k[1],simplex,theta,threads,progressbar)
  }
  y_xmap_x = RcppGCCM4Grid(effectmat,causemat,libsizes,pred,E[2],tau[2],k[2],simplex,theta,threads,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' geographical convergent cross mapping
#'
#' @param data The observation data.
#' @param cause Name of causal variable.
#' @param effect Name of effect variable.
#' @param libsizes A vector of library sizes to use.
#' @param E (optional) Dimensions of the embedding.
#' @param tau (optional) Step of spatial lags.
#' @param k (optional) Number of nearest neighbors to use for prediction.
#' @param theta (optional) Weighting parameter for distances, useful when `algorithm` is `smap`.
#' @param algorithm (optional) Algorithm used for prediction.
#' @param pred (optional) Row numbers(`vector` for lattice data) or row-column numbers(`matrix` for grid data) used for predictions.
#' @param nb (optional) The neighbours list.
#' @param threads (optional) Number of threads.
#' @param bidirectional (optional) whether to identify bidirectional causal associations.
#' @param trend.rm (optional) Whether to remove the linear trend.
#' @param progressbar (optional) whether to print the progress bar.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping prediction results}
#' \item{\code{varname}}{names of causal and effect variable}
#' \item{\code{bidirectional}}{whether to identify bidirectional causal associations}
#' }
#' @export
#' @importFrom methods setGeneric
#' @importFrom methods setMethod
#' @name gccm
#' @rdname gccm
#' @aliases gccm,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData"))
#' \donttest{
#' g = gccm(columbus,"HOVAL","CRIME",libsizes = seq(5,45,5),E = 6)
#' g
#' plot(g, ylimits = c(0,0.85))
#' }
methods::setMethod("gccm", "sf", .gccm_sf_method)

#' @rdname gccm
methods::setMethod("gccm", "SpatRaster", .gccm_spatraster_method)
