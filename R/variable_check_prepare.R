.check_character = \(x,...){
  xstrs = c(x,...)
  for (i in xstrs){
    if (!inherits(i,"character")) {
      stop("Please check the characters in the function input.")
    }
  }
  return(xstrs)
}

.check_inputelementnum = \(x,n,med.num = NULL){
  if (is.null(med.num) || length(x) == 1){
    res = rep(x,length.out = n)
  } else if (length(x) == 2) {
    res = c(rep(x[1],2),rep(x[2],med.num))
  } else {
    res = c(x[1:2],rep(x[c(-1,-2)],length.out = med.num))
  }
  return(abs(res))
}

.check_parallellevel = \(parallel.level){
  pl = 0
  if (parallel.level != "low"){
    pl = 1
  }
  return(pl)
}

.internal_varname = \(mediator = NULL){
  .varname = c("cause","effect")
  if (!is.null(mediator)){
    .varname = c(.varname,paste0("z",seq_along(mediator)))
  }
  return(.varname)
}

.internal_samplemat = \(mat,size = NULL,seed = 123){
  nnaindice = which(!is.na(mat), arr.ind = TRUE)
  if (is.null(size)){
    return(nnaindice)
  } else {
    set.seed(seed)
    indices = sample(nrow(nnaindice), size = min(size,nrow(nnaindice)), replace = FALSE)
    return(nnaindice[indices])
  }
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
