methods::setGeneric("gccm", function(data, ...) standardGeneric("gccm"))

.gccm_sf_method = \(data, cause, effect, libsizes, E = 3, tau = 1, k = E+1, theta = 1, algorithm = "simplex",
                    nb = NULL, threads = detectThreads(), trendRM = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  coords = sdsfun::sf_coordinates(data)
  cause = data[,cause,drop = TRUE]
  effect = data[,effect,drop = TRUE]
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  if (length(cause) != length(nb)) stop("Incompatible Data Dimensions!")
  if (trendRM){
    dtf = data.frame(cause = cause, effect = effect, x = coords[,1], y = coords[,2])
    cause = sdsfun::rm_lineartrend("cause~x+y", data = dtf)
    effect = sdsfun::rm_lineartrend("effect~x+y", data = dtf)
  }

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = RcppGCCM4Lattice(cause,effect,nb,libsizes,E,tau,k,simplex,theta,threads,progressbar)
  y_xmap_x = RcppGCCM4Lattice(effect,cause,nb,libsizes,E,tau,k,simplex,theta,threads,progressbar)

  return(.bind_xmapdf(x_xmap_y,y_xmap_x,varname))
}

.gccm_spatraster_method = \(data, cause, effect, libsizes, E = 3, tau = 1, k = E+3, theta = 1, algorithm = "simplex",
                            RowCol = NULL, threads = detectThreads(), trendRM = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  data = data[[c(cause,effect)]]
  names(data) = c("cause","effect")

  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (trendRM){
    dtf$cause = sdsfun::rm_lineartrend("cause~x+y", data = dtf)
    dtf$effect = sdsfun::rm_lineartrend("effect~x+y", data = dtf)
  }
  causemat = matrix(dtf[,3],nrow = terra::nrow(data),byrow = TRUE)
  effectmat = matrix(dtf[,4],nrow = terra::nrow(data),byrow = TRUE)

  maxlibsize = min(dim(causemat))
  selvec = seq(5,maxlibsize,5)
  if (is.null(RowCol)) RowCol = as.matrix(expand.grid(selvec,selvec))

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = RcppGCCM4Grid(causemat,effectmat,libsizes,RowCol,E,tau,k,simplex,theta,threads,progressbar)
  y_xmap_x = RcppGCCM4Grid(effectmat,causemat,libsizes,RowCol,E,tau,k,simplex,theta,threads,progressbar)

  return(.bind_xmapdf(x_xmap_y,y_xmap_x,varname))
}

#' geographical convergent cross mapping
#'
#' @param data The observation data.
#' @param cause Name of causal variable.
#' @param effect Name of effect variable.
#' @param libsizes A vector of library sizes to use.
#' @param E (optional) The dimensions of the embedding.
#' @param tau (optional) The step of spatial lags.
#' @param k (optional) Number of nearest neighbors to use for prediction.
#' @param theta (optional) Weighting parameter for distances, useful when `algorithm` is `smap`.
#' @param algorithm (optional) Algorithm used for prediction.
#' @param nb (optional) The neighbours list.
#' @param RowCol (optional) Matrix of selected row and cols numbers.
#' @param threads (optional) Number of threads.
#' @param trendRM (optional) Whether to remove the linear trend.
#' @param progressbar (optional) whether to print the progress bar.
#'
#' @return A list.
#' \describe{
#' \item{\code{xmap}}{cross-mapping prediction outputs}
#' \item{\code{varname}}{names of causal and effect variable}
#' }
#' @export
#' @importFrom methods setGeneric
#' @importFrom methods setMethod
#' @name gccm
#' @rdname gccm
#' @aliases gccm,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#'                        quiet=TRUE)
#' \donttest{
#' g = gccm(columbus, "HOVAL", "CRIME", libsizes = seq(5,45,5))
#' g
#' plot(g, ylimits = c(0,0.65))
#' }
methods::setMethod("gccm", "sf", .gccm_sf_method)

#' @rdname gccm
methods::setMethod("gccm", "SpatRaster", .gccm_spatraster_method)
