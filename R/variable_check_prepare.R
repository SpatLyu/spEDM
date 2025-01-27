.check_tgcharacter = \(tg){
  if (!inherits(tg,"character")) {
    stop("`target` must be character.")
  }
  return(tg)
}

.check_cecharacter = \(cause,effect){
  if (!inherits(cause,"character") || !inherits(effect,"character")) {
    stop("The `cause` and `effect` must be character.")
  }
  varname = c(cause,effect)
  return(varname)
}

.check_input2element = \(x){
  if (length(x) == 1) x = rep(x,times = 2)
  if (length(x) > 2) x = x[c(1,2)]
  return(x)
}

.check_indices = \(x,totalnum){
  return(x[(x<=totalnum)&(x>=1)])
}

.uni_lattice = \(data,target){
  target = .check_tgcharacter(target)
  res = data[,target,drop = TRUE]
  return(res)
}

.uni_grid = \(data,target){
  target = .check_tgcharacter(target)
  data = data[[target]]
  res = matrix(terra::values(data),nrow = terra::nrow(data),byrow = TRUE)
  return(res)
}
