#' Internal: is `x` a list of per-control-cell objects?
#'
#' Several functions accept either a single object (tabulation, Pi, etc.) or
#' a list of such objects split by control cells. The old test compared
#' `class(x)` to the string "list", which breaks for multi-class objects (e.g. data.tables,
#' where the length-2 class vector makes `if()` error) and misfires for
#' data.frames, which are lists. This helper captures the intended check.
#'
#' @noRd
is_cell_list = function(x){
  is.list(x) && !is.data.frame(x)
}
