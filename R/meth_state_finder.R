#' Mixture Coefficient Estimation
#'
#' This function determines if a CpG site is differentially methylated using
#' Hidden Markov Models.
#' States:
#' hypo-methylated  = 1;
#' non-differential = 2;
#' hyper-methylated = 3.
#'
#' @param input_signal A numeric vector of AUC scores
#' @param input_pos An integer vector of chromosomal locations
#' @param auc_sd Standard deviation of AUC signal
#' @param PT_start Fraction of AUC signal lower than 0.1 or greater than 0.9
#' @param length_cutoff An integer to remove streches of differential
#' methylation shorter than cutoff
#' @return An integer vector of methylated states
#' @export
meth_state_finder <- function(input_signal, input_pos, auc_sd, PT_start,
  length_cutoff) {
  # check parameters
  stopifnot(all(!is.na(input_signal)))
  stopifnot(is.numeric(input_signal))
  input_pos <- as.integer(input_pos)
  length_cutoff <- as.integer(length_cutoff)

  # omit NAs
  final_states <- rep(NA, length(input_signal))
  idx_not_na <- which(!is.na(input_signal))
  input_signal <- input_signal[idx_not_na]

  muk <- c(0.25, 0.5, 0.75)
  NormDist <- 1e5
  ratioSD <- 0.4
  sepsilon <- rep(auc_sd * ratioSD, length(muk))
  sepsilon[which(muk == 0.5)] <- auc_sd * (1 - ratioSD)

  KS <- length(muk)
  CovPos <- diff(input_pos)
  CovDist <- CovPos/NormDist
  CovDist1 <- log(1 - exp(-CovDist))
  W <- length(input_signal)
  NCov <- length(CovDist)
  PT <- log(rep(PT_start, KS))

  P <- matrix(data = 0, nrow = KS, ncol = (KS * NCov))
  emission <- matrix(data = 0, nrow = KS, ncol = W)

  #### Calculates Transition and Emission Probabilities ##
  out <- .Fortran("transemisi", as.vector(muk), as.integer(NCov),
    as.vector(input_signal), as.integer(KS), as.vector(CovDist1),
    as.vector(sepsilon), as.integer(W), as.matrix(PT), as.matrix(P),
    as.matrix(emission))
  P <- out[[9]]
  emission <- out[[10]]
  etav <- log(rep(1, KS) * (1/KS))
  psi <- matrix(data = 0, nrow = KS, ncol = W)
  path <- c(as.integer(rep(0, W)))

  ##### Viterbi Algorithm ####
  out2 <- .Fortran("bioviterbii", as.vector(etav), as.matrix(P),
    as.matrix(emission), as.integer(W), as.integer(KS), as.vector(path),
    as.matrix(psi))
  meth_states <- out2[[6]]

  # fix output: "undifferentiate" short segments, reinsert NAs, correct NAs
  final_states <- fix_short_segments(meth_states, cutoff = length_cutoff)
  return(final_states)
}

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
