#' detect the max number of threads that can be used
#'
#' @return An integer
#' @export
#'
#' @examples
#' \donttest{
#' detectThreads()
#' }
detectThreads = \() {
  return(DetectMaxNumThreads())
}
