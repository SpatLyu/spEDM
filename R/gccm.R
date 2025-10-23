.gccm_sf_method = \(data, cause, effect, libsizes = NULL, E = 3, tau = 1, k = E+2, theta = 1, algorithm = "simplex", lib = NULL, pred = NULL, style = 1, stack = FALSE,
                    dist.metric = "L2",dist.average = TRUE,nb = NULL,threads = detectThreads(),parallel.level = "low",bidirectional = TRUE,detrend = TRUE,progressbar = TRUE){
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

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCCM4Lattice(cause,effect,nb,libsizes,lib,pred,E[1],tau[1],k[1],simplex,theta,threads,pl,
                                style, stack, .check_distmetric(dist.metric), dist.average, TRUE, progressbar)
  }
  y_xmap_x = RcppGCCM4Lattice(effect,cause,nb,libsizes,lib,pred,E[2],tau[2],k[2],simplex,theta,threads,pl,
                              style, stack, .check_distmetric(dist.metric), dist.average, TRUE, progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

.gccm_spatraster_method = \(data, cause, effect, libsizes = NULL, E = 3, tau = 1, k = E+2, theta = 1, algorithm = "simplex", lib = NULL, pred = NULL, style = 1, stack = FALSE, dist.metric = "L2",
                            dist.average = TRUE, embed.direction = 0, win.ratio = 0, threads = detectThreads(), parallel.level = "low", bidirectional = TRUE, detrend = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  tau = .check_inputelementnum(tau,2)
  k = .check_inputelementnum(k,2)
  win.ratio = .check_inputelementnum(win.ratio,2)
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

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCCM4Grid(causemat,effectmat,libsizes,lib,pred,E[1],tau[1],k[1],simplex,theta,threads,pl,style,stack,
                             .check_distmetric(dist.metric),dist.average,TRUE,embed.direction,win.ratio,progressbar)
  }
  y_xmap_x = RcppGCCM4Grid(effectmat,causemat,libsizes,lib,pred,E[2],tau[2],k[2],simplex,theta,threads,pl,style,stack,
                           .check_distmetric(dist.metric),dist.average,TRUE,embed.direction,win.ratio,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' geographical convergent cross mapping
#'
#' @inheritParams gcmc
#' @param theta (optional) weighting parameter for distances, useful when `algorithm` is `smap`.
#' @param algorithm (optional) prediction algorithm.
#' @param stack (optional) whether to stack embeddings.
#' @param dist.average (optional) whether to average distance.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping results}
#' \item{\code{varname}}{names of causal and effect variable}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name gccm
#' @aliases gccm,sf-method
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
#' @param embed.direction (optional) direction selector for embeddings (`0` returns all directions, `1-8` correspond to NW, N, NE, W, E, SW, S, SE).
#' @param win.ratio (optional) ratio of sliding window scale to speed up state-space predictions.
methods::setMethod("gccm", "SpatRaster", .gccm_spatraster_method)
