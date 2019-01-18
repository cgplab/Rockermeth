#' Fix Short Methylation States
#'
#' Given a vector of methylation states, this function sets any hyper-(3) or
#' hypo-(1) methylated stretches of sites shorter than or equal to cutoff, to
#' non-differential state (2).
#'
#' @param meth_states An integer vector of methylation states
#' @param cutoff Length cutoff: longer segments will be ignored
#' @keywords internal
fix_short_segments <- function(meth_states, cutoff) {
  assertthat::assert_that(all(meth_states %in% 1:3))

  rle_states <- S4Vectors::Rle(meth_states)
  rle_length <- S4Vectors::runLength(rle_states)
  rle_value <- S4Vectors::runValue(rle_states)
  starts <- S4Vectors::start(rle_states)
  ends <- S4Vectors::end(rle_states)

  idx <- which(rle_length <= cutoff & rle_value != 2)
  for (i in idx) {
    meth_states[starts[i]:ends[i]] <- 2
  }
  return(meth_states)
}
