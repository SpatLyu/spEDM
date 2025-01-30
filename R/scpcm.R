.scpcm_sf_method = \(data, cause, effect, libsizes, E = c(3,3), tau = 1, k = 4, theta = 1, algorithm = "simplex",
                     nb = NULL, threads = detectThreads(), include.self = FALSE, trendRM = TRUE, progressbar = TRUE){
  varname = .check_cecharacter(cause, effect)
  E = .check_input2element(E)
  k = .check_input2element(k)
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
  x_xmap_y = RcppGCCM4Lattice(cause,effect,nb,libsizes,E[1],tau,k[1],simplex,theta,threads,include.self,progressbar)
  y_xmap_x = RcppGCCM4Lattice(effect,cause,nb,libsizes,E[2],tau,k[2],simplex,theta,threads,include.self,progressbar)

  return(.bind_xmapdf(x_xmap_y,y_xmap_x,varname))
}