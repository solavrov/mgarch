
#' return sublist of a list
#'
#' @param l given list
#' @param n1 start index of sublist
#' @param n2 end index of sublist
#'
#' @return sublist
#' @export
#'
#' @examples
helper.sublist <- function(l, n1, n2=length(l)) {
  l2 <- list()
  for (i in n1:n2) l2 <- append(l2, l[[i]])
  return (l2)
}
