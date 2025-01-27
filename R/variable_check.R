.check_character = \(cause,effect){
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
