.check_character = \(x,...){
  xstrs = c(x,...)
  for (i in xstrs){
    if (!inherits(i,"character")) {
      stop("Please check the characters in the function input.")
    }
  }
  return(xstrs)
}

.check_inputelementnum = \(x,n){
  return(abs(rep(x,length.out = n)))
}

.check_indices = \(x,totalnum){
  return(x[(x<=totalnum)&(x>=1)])
}

.check_libsizes_gird = \(libsizes){
  if (!inherits(libsizes,"matrix")){
    libsizes = matrix(rep(libsizes,times = 2),ncol = 2)
  }
  return(libsizes)
}

.internal_varname = \(mediator = NULL){
  .varname = c("cause","effect")
  if (!is.null(mediator)){
    .varname = c(.varname,paste0("z",seq_along(mediator)))
  }
  return(.varname)
}

.internal_predmat = \(mat){
  maxlibsize = min(dim(mat))
  selvec = seq(5,maxlibsize,5)
  pred = as.matrix(expand.grid(selvec,selvec))
  return(pred)
}

.internal_lattice_nb = \(data){
  if (sdsfun::sf_geometry_type(data) %in% c('point','multipoint')){
    nb = sdsfun::spdep_nb(data,k = 8)
  } else {
    nb = sdsfun::spdep_nb(data)
  }
  return(nb)
}

.internal_trend_rm = \(data,.varname,coords = NULL){
  if (is.null(coords)){
    for (i in seq_along(.varname)){
      data[,.varname[i]] = sdsfun::rm_lineartrend(paste0(.varname[i],"~x+y"), data = data)
    }
  } else {
    data = dplyr::bind_cols(data,coords)
    for (i in seq_along(.varname)){
      data[,.varname[i]] = sdsfun::rm_lineartrend(paste0(.varname[i],"~X+Y"), data = data)
    }
  }
  return(data)
}

.uni_lattice = \(data,target,trend.rm){
  target = .check_character(target)
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  data = sf::st_drop_geometry(data)
  data = data[,target,drop = FALSE]
  names(data) = "target"
  if (trend.rm){
    data = .internal_trend_rm(data,"target",coords)
  }
  res = data[,"target",drop = TRUE]
  return(res)
}

.uni_grid = \(data,target,trend.rm){
  target = .check_character(target)
  data = data[[target]]
  names(data) = "target"
  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (trend.rm){
    dtf = .internal_trend_rm(dtf,"target")
  }
  res = matrix(dtf[,"target"],nrow = terra::nrow(data),byrow = TRUE)
  return(res)
}

.multivar_lattice = \(data,columns,trend.rm){
  columns = .check_character(columns)
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  data = sf::st_drop_geometry(data)
  data = data[,columns,drop = FALSE]
  .varname = paste0("z",seq_along(columns))
  names(data) = .varname
  if (trend.rm){
    data = .internal_trend_rm(data,.varname,coords)
  }
  res = as.matrix(data[,.varname,drop = FALSE])
  return(res)
}

.multivar_grid = \(data,columns,trend.rm){
  columns = .check_character(columns)
  data = data[[columns]]
  .varname = paste0("z",seq_along(columns))
  names(data) = .varname
  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (trend.rm){
    dtf = .internal_trend_rm(dtf,.varname)
  }
  res = as.matrix(dtf[,.varname,drop = FALSE])
  return(res)
}
