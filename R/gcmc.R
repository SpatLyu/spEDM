methods::setGeneric("gcmc", function(data, ...) standardGeneric("gcmc"))

.gcmc_sf_method = \(data, cause, effect, E = c(3,3), k = 4, r = k + 10, nb = NULL,
                    threads = detectThreads(), include.self = FALSE, trendRM = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  k = .check_inputelementnum(k,2)
  r = .check_inputelementnum(r,2)
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

  res = RcppGCMC4Lattice(cause,effect,nb,E,k,r,threads,include.self,progressbar)
  names(res) = .name_xmap2cause(varname)
  return(res)
}

.gcmc_spatraster_method = \(data, cause, effect, E = c(3,3), k = 4, r = k + 10, threads = detectThreads(),
                            include.self = FALSE, trendRM = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  k = .check_inputelementnum(k,2)
  r = .check_inputelementnum(r,2)
  data = data[[c(cause,effect)]]
  names(data) = c("cause","effect")

  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (trendRM){
    dtf$cause = sdsfun::rm_lineartrend("cause~x+y", data = dtf)
    dtf$effect = sdsfun::rm_lineartrend("effect~x+y", data = dtf)
  }
  causemat = matrix(dtf[,"cause"],nrow = terra::nrow(data),byrow = TRUE)
  effectmat = matrix(dtf[,"effect"],nrow = terra::nrow(data),byrow = TRUE)

  res = RcppGCMC4Grid(causemat,effectmat,E,k,r,threads,include.self,progressbar)
  names(res) = .name_xmap2cause(varname)
  return(res)
}

#' geographical cross mapping cardinality
#'
#' @param data The observation data.
#' @param cause Name of causal variable.
#' @param effect Name of effect variable.
#' @param E (optional) The dimensions of the embedding.
#' @param k (optional) Number of nearest neighbors to use for prediction.
#' @param r (optional) Maximum number of neighbors usable for intersection cardinality computation.
#' @param algorithm (optional) Algorithm used for prediction.
#' @param nb (optional) The neighbours list.
#' @param threads (optional) Number of threads.
#' @param include.self (optional) Whether to include the current state when constructing the embedding vector.
#' @param trendRM (optional) Whether to remove the linear trend.
#' @param progressbar (optional) whether to print the progress bar.
#'
#' @return A numeric vector.
#'
#' @export
#' @importFrom methods setGeneric
#' @importFrom methods setMethod
#' @name gcmc
#' @rdname gcmc
#' @aliases gcmc,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#'                        quiet=TRUE)
#' \donttest{
#' g = gcmc(columbus,"HOVAL","CRIME",E = c(6,5))
#' g
#' }
methods::setMethod("gcmc", "sf", .gcmc_sf_method)

#' @rdname gcmc
methods::setMethod("gcmc", "SpatRaster", .gcmc_spatraster_method)
