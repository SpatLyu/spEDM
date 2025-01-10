.check_character = \(cause,effect){
  if (!inherits(cause,"character") || !inherits(effect,"character")) {
    stop("The `cause` and `effect` must be character.")
  }
  varname = c(cause,effect)
  return(varname)
}
