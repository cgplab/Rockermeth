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
#' @param length_cutoff An integer to remove streches of differential
#' methylation shorter than cutoff
#' @param na_cutoff An integer to set stretches of NAs to the state of flanking sites
#' depending if the stretches are shorter than cutoff
#'
#' @return An integer vector of methylated states
#' @export

meth_state_finder <- function(input_signal, input_pos, length_cutoff, na_cutoff) {
  stopifnot(is.numeric(input_signal))
  input_pos <- as.integer(input_pos)
  length_cutoff <- as.integer(length_cutoff)
  na_cutoff <- as.integer(na_cutoff)

  final_states <- rep(NA, length(input_signal))
  idx_not_na <- which(!is.na(input_signal))
  input_signal <- input_signal[idx_not_na]

  muk <- c(0.25, 0.5, 0.75)
  NormDist <- 1e+05
  PTStart <- 0.05
  ratioSD <- 0.4
  SDtot <- sd(input_signal)
  sepsilon <- rep(SDtot * ratioSD, length(muk))
  sepsilon[which(muk == 0.5)] <- SDtot * (1 - ratioSD)

  KS <- length(muk)
  CovPos <- diff(input_pos)
  CovDist <- CovPos/NormDist
  CovDist1 <- log(1 - exp(-CovDist))
  W <- length(input_signal)
  NCov <- length(CovDist)
  PT <- log(rep(PTStart, KS))
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

  # fix output: "undifferentiate" short segments, restore NAs, correct NAs
  final_states[idx_not_na] <- fix_short_segments(meth_states, cutoff = length_cutoff)
  final_states <- fix_na_segments(final_states, cutoff = na_cutoff)
  return(final_states)
}

#' Fix Short Methylation States
#'
#' Given a vector of methylation states, this function sets to non-differential
#' (2) any hyper-(3) or hypo-(1) methylated stretches of sites shorter than
#' or equal to cutoff.
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

#' Fix NAs Methylation States
#'
#' Given a vector of methylation states, this function sets NA stretches
#' shorter than ore equal to cutoff to the state of flanking sites if flanking
#' sites have the same state.
#'
#' @inheritParams fix_short_segments
#' @keywords internal
fix_na_segments <- function(meth_states, cutoff) {
  rle_states <- S4Vectors::Rle(meth_states)
  rle_length <- S4Vectors::runLength(rle_states)
  rle_value <- S4Vectors::runValue(rle_states)
  starts <- S4Vectors::start(rle_states)
  ends <- S4Vectors::end(rle_states)
  N <- S4Vectors::nrun(rle_states)

  idx_na <- is.na(rle_value)
  idx_short <- rle_length <= cutoff
  idx_same_flanks <- c(F, rle_value[1:(N-2)] == rle_value[3:N], F)
  idx <- which(idx_na & idx_short & idx_same_flanks)
  # set NA state to state of previous stretch (equal to state of the following one)
  for (i in idx){
    meth_states[starts[i]:ends[i]] <- rle_value[i - 1]
  }
  return(meth_states)
}
