#' methods::setGeneric("scpcm", function(data, ...) standardGeneric("scpcm"))
#'
#' .scpcm_sf_method = \(data, cause, effect, mediator, libsizes, E = 3, tau = 1, k = 4, theta = 1, algorithm = "simplex", nb = NULL,
#'                      threads = detectThreads(), bidirectional = TRUE, cumulate = FALSE, include.self = TRUE, trendRM = TRUE, progressbar = TRUE){
#'   varname = .check_character(c(cause, effect, mediator))
#'   E = .check_inputelementnum(E,length(varname))
#'   k = .check_inputelementnum(k,2)
#'   coords = as.data.frame(sdsfun::sf_coordinates(data))
#'   if (is.null(nb)) nb = sdsfun::spdep_nb(data)
#'   if (nrow(data) != length(nb)) stop("Incompatible Data Dimensions!")
#'   data = sf::st_drop_geometry(data)
#'   data = data[,varname]
#'
#'   if (trendRM){
#'     data = dplyr::bind_cols(data,coords)
#'     for (i in seq_along(varname)){
#'       data[,varname[i]] = sdsfun::rm_lineartrend(paste0(varname[i],"~X+Y"), data = data)
#'     }
#'   }
#'
#'   cause = data[,cause,drop = TRUE]
#'   effect = data[,effect,drop = TRUE]
#'   medmat = as.matrix(data[,mediator])
#'
#'   simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
#'   x_xmap_y = NULL
#'   if (bidirectional){
#'     x_xmap_y = RcppSCPCM4Lattice(cause,effect,medmat,nb,libsizes,E[-2],tau,k[1],simplex,theta,threads,cumulate,include.self,progressbar)
#'   }
#'   y_xmap_x = RcppSCPCM4Lattice(effect,cause,medmat,nb,libsizes,E[-1],tau,k[2],simplex,theta,threads,cumulate,include.self,progressbar)
#'
#'   return(.bind_xmapdf2(varname,x_xmap_y,y_xmap_x,bidirectional))
#' }
#'
#' .scpcm_spatraster_method = \(data, cause, effect, mediator, libsizes, E = 3, tau = 1, k = 4, theta = 1, algorithm = "simplex", RowCol = NULL,
#'                              threads = detectThreads(), bidirectional = TRUE, cumulate = FALSE, include.self = TRUE, trendRM = TRUE, progressbar = TRUE){
#'   varname = .check_character(cause, effect, mediator)
#'   E = .check_inputelementnum(E,length(varname))
#'   k = .check_inputelementnum(k,2)
#'   data = data[[varname]]
#'   zs = paste0("z",seq_along(mediator))
#'   names(data) = c(c("cause","effect"),zs)
#'
#'   dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
#'   if (trendRM){
#'     for (i in seq_along(varname)){
#'       dtf[,varname[i]] = sdsfun::rm_lineartrend(paste0(varname[i],"~x+y"), data = dtf)
#'     }
#'   }
#'   causemat = matrix(dtf[,"cause"],nrow = terra::nrow(data),byrow = TRUE)
#'   effectmat = matrix(dtf[,"effect"],nrow = terra::nrow(data),byrow = TRUE)
#'   medmat = as.matrix(dtf[,zs])
#'
#'   maxlibsize = min(dim(causemat))
#'   selvec = seq(5,maxlibsize,5)
#'   if (is.null(RowCol)) RowCol = as.matrix(expand.grid(selvec,selvec))
#'
#'   simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
#'   x_xmap_y = NULL
#'   if (bidirectional){
#'     x_xmap_y = RcppSCPCM4Grid(causemat,effectmat,medmat,libsizes,E[-2],RowCol,tau,k[1],simplex,theta,threads,cumulate,include.self,progressbar)
#'   }
#'   y_xmap_x = RcppSCPCM4Grid(effectmat,causemat,medmat,libsizes,E[-1],RowCol,tau,k[2],simplex,theta,threads,cumulate,include.self,progressbar)
#'
#'   return(.bind_xmapdf2(varname,x_xmap_y,y_xmap_x,bidirectional))
#' }
#'
#' #' spatially convergent partial cross mapping
#' #'
#' #' @inheritParams gccm
#' #' @param mediator Name of mediator variable.
#' #' @param cumulate (optional) Serial or cumulative computation of partial cross mapping.
#' #'
#' #' @return A list.
#' #' \describe{
#' #' \item{\code{pxmap}}{partial cross mapping prediction results}
#' #' \item{\code{xmap}}{cross mapping prediction results}
#' #' \item{\code{varname}}{names of causal and effect variable}
#' #' \item{\code{bidirectional}}{whether to identify bidirectional potential causal relationships}
#' #' }
#' #' @name scpcm
#' #' @rdname scpcm
#' #' @aliases scpcm,sf-method
#' #'
#' #' @examples
#' #' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#' #'                        quiet=TRUE)
#' #' \dontrun{
#' #' g = scpcm(columbus, "HOVAL", "CRIME", "INC",
#' #'           libsizes = seq(5,40,5), E = c(6,5,3))
#' #' g
#' #' plot(g, ylimits = c(0,0.8))
#' #' }
#' methods::setMethod("scpcm", "sf", .scpcm_sf_method)
#'
#' #' @rdname scpcm
#' methods::setMethod("scpcm", "SpatRaster", .scpcm_spatraster_method)
