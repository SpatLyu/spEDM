methods::setGeneric("gccm", function(data, ...) standardGeneric("gccm"))

.gccm_sf_method = \(data, cause, effect, libsizes, E = c(3,3), tau = 1, k = 4, theta = 1, algorithm = "simplex", nb = NULL,
                    threads = detectThreads(), bidirectional = TRUE, include.self = FALSE, trendRM = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  k = .check_inputelementnum(k,2)
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  if (nrow(data) != length(nb)) stop("Incompatible Data Dimensions!")
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  data = sf::st_drop_geometry(data)
  data = data[,varname]

  if (trendRM){
    data = dplyr::bind_cols(data,coords)
    for (i in 1:length(varname)){
      data[,varname[i]] = sdsfun::rm_lineartrend(paste0(varname[i],"~X+Y"), data = data)
    }
  }

  cause = data[,cause,drop = TRUE]
  effect = data[,effect,drop = TRUE]

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCCM4Lattice(cause,effect,nb,libsizes,E[1],tau,k[1],simplex,theta,threads,include.self,progressbar)
  }
  y_xmap_x = RcppGCCM4Lattice(effect,cause,nb,libsizes,E[2],tau,k[2],simplex,theta,threads,include.self,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

.gccm_spatraster_method = \(data, cause, effect, libsizes, E = c(3,3), tau = 1, k = 4, theta = 1, algorithm = "simplex", RowCol = NULL,
                            threads = detectThreads(), bidirectional = TRUE, include.self = FALSE, trendRM = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  k = .check_inputelementnum(k,2)
  data = data[[c(cause,effect)]]
  names(data) = c("cause","effect")

  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (trendRM){
    dtf$cause = sdsfun::rm_lineartrend("cause~x+y", data = dtf)
    dtf$effect = sdsfun::rm_lineartrend("effect~x+y", data = dtf)
  }
  causemat = matrix(dtf[,"cause"],nrow = terra::nrow(data),byrow = TRUE)
  effectmat = matrix(dtf[,"effect"],nrow = terra::nrow(data),byrow = TRUE)

  maxlibsize = min(dim(causemat))
  selvec = seq(5,maxlibsize,5)
  if (is.null(RowCol)) RowCol = as.matrix(expand.grid(selvec,selvec))

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCCM4Grid(causemat,effectmat,libsizes,RowCol,E[1],tau,k[1],simplex,theta,threads,include.self,progressbar)
  }
  y_xmap_x = RcppGCCM4Grid(effectmat,causemat,libsizes,RowCol,E[2],tau,k[2],simplex,theta,threads,include.self,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
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
#' @param bidirectional (optional) whether to identify bidirectional potential causal relationships.
#' @param include.self (optional) Whether to include the current state when constructing the embedding vector.
#' @param trendRM (optional) Whether to remove the linear trend.
#' @param progressbar (optional) whether to print the progress bar.
#'
#' @return A list.
#' \describe{
#' \item{\code{xmap}}{cross mapping prediction outputs}
#' \item{\code{varname}}{names of causal and effect variable}
#' \item{\code{bidirectional}}{whether to identify bidirectional potential causal relationships}
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
#' g = gccm(columbus,"HOVAL","CRIME",libsizes = seq(5,40,5),E = c(6,5))
#' g
#' plot(g, ylimits = c(0,0.8))
#' }
methods::setMethod("gccm", "sf", .gccm_sf_method)

#' @rdname gccm
methods::setMethod("gccm", "SpatRaster", .gccm_spatraster_method)
