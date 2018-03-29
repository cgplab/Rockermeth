#' Mixture Coefficient Estimation
#'
#' This function classifies CpG sites into hypo-methylated, not differentially
#' methylated and hyper-methylated states using a Heterogeneous Shifting Level
#' Model (HSLM).
#' States:
#' hypo-methylated  = 1;
#' non-differential = 2;
#' hyper-methylated = 3.
#'
#' @param input_signal A numeric vector of AUC scores
#' @param input_pos An integer vector of chromosomal locations
#' @param auc_sd Standard deviation of AUC signal (genome wide)
#' @param length_cutoff An integer to remove streches of differential
#' methylation shorter than cutoff
#' @param pt_start Transition probability of the HSLM. Default is 0.05.
#' @param normdist Distance normalization parameter of the HSLM. Default is 1e5.
#' @param ratiosd Fraction between the standard deviation of AUC values of
#' differentially methylated sites and the total standard deviation of AUC
#' values. Default is 0.4.
#' @return An integer vector of the states of methylation.
#' @export
meth_state_finder <- function(input_signal, input_pos, auc_sd, length_cutoff,
                              pt_start, normdist, ratiosd) {
  assertthat::assert_that(all(!is.na(input_signal)))
  assertthat::assert_that(is.numeric(input_signal))
  assertthat::assert_that(length(input_signal) == length(input_pos))

  input_pos <- as.integer(input_pos)
  length_cutoff <- as.integer(length_cutoff)

  muk <- c(0.25, 0.5, 0.75)
  sepsilon <- rep(auc_sd * ratiosd, length(muk))
  sepsilon[which(muk == 0.5)] <- auc_sd * (1 - ratiosd)

  KS <- length(muk)
  CovPos <- diff(input_pos)
  CovDist <- CovPos/normdist
  CovDist1 <- log(1 - exp(-CovDist))
  W <- length(input_signal)
  NCov <- length(CovDist)
  PT <- log(rep(pt_start, KS))

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
  path <- as.integer(rep(0, W))

  ##### Viterbi Algorithm ####
  out2 <- .Fortran("bioviterbii", as.vector(etav), as.matrix(P),
                   as.matrix(emission), as.integer(W), as.integer(KS), as.vector(path),
                   as.matrix(psi))
  meth_states <- out2[[6]]

  # fix output: "undifferentiate" short segments, reinsert NAs, correct NAs
  return(fix_short_segments(meth_states, cutoff = length_cutoff))
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
  assertthat::assert_that(all(meth_states >= 1 & meth_states <= 3))
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
