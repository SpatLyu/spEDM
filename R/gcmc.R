methods::setGeneric("gcmc", function(data, ...) standardGeneric("gcmc"))

.gcmc_sf_method = \(data, cause, effect, E = c(3,3), tau = 1, k = 4, r = k + 10, pred = NULL,
                    nb = NULL, threads = detectThreads(), trend.rm = FALSE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  k = .check_inputelementnum(k,2)
  r = .check_inputelementnum(r,2)
  tau = .check_inputelementnum(tau,2)
  .varname = .internal_varname()
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  if (nrow(data) != length(nb)) stop("Incompatible Data Dimensions!")
  if (is.null(pred)) pred = 1:nrow(data)
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  data = sf::st_drop_geometry(data)
  data = data[,varname]
  names(data) = .varname

  if (trend.rm){
    data = dplyr::bind_cols(data,coords)
    for (i in seq_along(.varname)){
      data[,.varname[i]] = sdsfun::rm_lineartrend(paste0(.varname[i],"~X+Y"), data = data)
    }
  }

  cause = data[,"cause",drop = TRUE]
  effect = data[,"effect",drop = TRUE]

  res = RcppGCMC4Lattice(cause,effect,nb,E,k,r,threads,include.self,progressbar)
  names(res) = .name_xmap2cause(varname)
  return(res)
}

.gcmc_spatraster_method = \(data, cause, effect, E = c(3,3), tau = 1, k = 4, r = k + 10, pred = NULL,
                            threads = detectThreads(),trend.rm = FALSE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  k = .check_inputelementnum(k,2)
  r = .check_inputelementnum(r,2)
  tau = .check_inputelementnum(tau,2)
  .varname = .internal_varname()
  data = data[[varname]]
  names(data) = .varname

  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (trend.rm){
    for (i in seq_along(.varname)){
      dtf[,.varname[i]] = sdsfun::rm_lineartrend(paste0(.varname[i],"~x+y"), data = dtf)
    }
  }
  causemat = matrix(dtf[,"cause"],nrow = terra::nrow(data),byrow = TRUE)
  effectmat = matrix(dtf[,"effect"],nrow = terra::nrow(data),byrow = TRUE)

  if (is.null(pred)) pred = .internal_predmat(causemat)

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
#' @param pred (optional) The row numbers(`vector`) of lattice data or the row-column numbers(`matrix`) of grid data used for predictions.
#' @param nb (optional) The neighbours list.
#' @param threads (optional) Number of threads.
#' @param trend.rm (optional) Whether to remove the linear trend.
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
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData"))
#' \donttest{
#' g = gcmc(columbus,"HOVAL","CRIME",E = c(6,5))
#' g
#' }
methods::setMethod("gcmc", "sf", .gcmc_sf_method)

#' @rdname gcmc
methods::setMethod("gcmc", "SpatRaster", .gcmc_spatraster_method)
