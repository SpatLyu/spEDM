#' detect the number of available threads
#'
#' @return An integer
#' @export
#'
#' @examples
#' detectThreads()
#' 
detectThreads = \() {
  return(DetectMaxNumThreads())
}
