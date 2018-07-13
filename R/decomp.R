#' @useDynLib rSeqSliceDecomp slicify
#' @export
decomp <- function(x)  {
  .Call(slicify, as.integer(x), as.integer(length(x))) 
}
