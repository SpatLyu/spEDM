.scpcm_sf_method = \(data, cause, effect, mediator, libsizes, E = 3, tau = 1, k = 4, theta = 1, algorithm = "simplex", nb = NULL,
                     threads = detectThreads(), bidirectional = TRUE, include.self = FALSE, trendRM = TRUE, progressbar = TRUE){
  varname = .check_character(c(cause, effect, mediator))
  E = .check_inputelementnum(E,length(varname))
  k = .check_inputelementnum(k,2)
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  if (nrow(data) != length(nb)) stop("Incompatible Data Dimensions!")
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
  medmat = as.matrix(data[,mediator])

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppSCPCM4Lattice(cause,effect,medmat,nb,libsizes,E[-2],tau,k[1],simplex,theta,threads,include.self,progressbar)
  }
  y_xmap_x = RcppSCPCM4Lattice(effect,cause,medmat,nb,libsizes,E[-1],tau,k[2],simplex,theta,threads,include.self,progressbar)

  return(.bind_xmapdf(varname[1:2],x_xmap_y,y_xmap_x,bidirectional))
}
