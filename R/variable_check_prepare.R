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
  if (length(x) == 1) x = rep(x,times = n)
  if (length(x) > 2) x = x[seq(1,n,by = 1)]
  return(x)
}

.check_indices = \(x,totalnum){
  return(x[(x<=totalnum)&(x>=1)])
}

.uni_lattice = \(data,target){
  target = .check_character(target)
  res = data[,target,drop = TRUE]
  return(res)
}

.uni_grid = \(data,target){
  target = .check_character(target)
  data = data[[target]]
  res = matrix(terra::values(data),nrow = terra::nrow(data),byrow = TRUE)
  return(res)
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
