#' Mixture Coefficient Estimation
#'
#' This function determines if a CpG site is differentially methylated using Hidden Markov Models.
#' States:
#' hypo-methylated = 1;
#' non-differential = 2;
#' hyper-methylated = 3.
#'
#' @param input_signal A numeric vector of AUC scores
#' @param input_pos An integer vector of chromosomal location
#' @param length_cutoff An integer to remove streches of methylation that are
#' too short
#' @param na_cutoff An integer to set NAs to flanking states depending on their length
#'
#' @return An integer vector of methylated states
#'
#' @export
#' @useDynLib MethLibrary

meth_state_finder <- function(input_signal, input_pos, length_cutoff = 5, na_cutoff = 2) {
  stopifnot(is.numeric(input_signal))
  input_pos <- as.integer(input_pos)
  length_cutoff <- as.integer(length_cutoff)
  na_cutoff <- as.integer(na_cutoff)

  final_states <- rep(NA, length(input_signal))
  idx_not_na <- which(!is.na(input_signal))
  input_signal <- input_signal[!is.na(input_signal)]

  muk <-  c(0.25, 0.5, 0.75)
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
  out <- .Fortran("transemisi", as.vector(muk), as.integer(NCov), as.vector(input_signal),
                  as.integer(KS), as.vector(CovDist1), as.vector(sepsilon), as.integer(W),
                  as.matrix(PT), as.matrix(P), as.matrix(emission))
  P <- out[[9]]
  emission <- out[[10]]
  etav <- log(rep(1, KS) * (1/KS))
  psi <- matrix(data = 0, nrow = KS, ncol = W)
  path <- c(as.integer(rep(0, W)))

  ##### Viterbi Algorithm ####
  out2 <- .Fortran("bioviterbii", as.vector(etav), as.matrix(P), as.matrix(emission),
                   as.integer(W), as.integer(KS), as.vector(path), as.matrix(psi))
  meth_states <- out2[[6]]

  # fix output
  final_states[idx_not_na] <- fix_short_segments(meth_states, cutoff = length_cutoff)
  final_states <- fix_na_states(final_states, cutoff = na_cutoff)
  return(final_states)
}

#' Fix Short Methylation States
#'
#' Given a vector of methylation states, this function sets to non-differential
#' (2) any hyper- (3) or hypo- (1) methylated stretches of sites with a length
#' lower than or equal to cutoff.
#'
#' @param meth_states An integer vector of methylation states
#' @param cutoff An integer

#' @keywords internal
fix_short_segments <- function(meth_states, cutoff) {
  rle_states <- S4Vectors::Rle(meth_states)
  idx <- which(S4Vectors::runLength(rle_states) < cutoff & S4Vectors::runValue(rle_states) != 2)
  for (i in idx) {
    meth_states[S4Vectors::start(rle_states)[i]:S4Vectors::end(rle_states)[i]] <- 2
  }
  return(meth_states)
}

#' Fix NAs Methylation States
#'
#' Change NA Methylation States according to flanking sites.
#' If a stretch of NAs has flanking sites with the same state,
#' the function changes the NAs to the state of the flaking sites. 
#' If flanking sites are different, if the stretch of NAs is shorter that cutoff, 
#' the function changes the NAs to state of the previous stretch (if first stretch, to the
#' state of following strech). 
#'
#' @param meth_states An integer vector of methylation states
#' @param cutoff An integer
#'
#' @keywords internal
fix_na_states <- function(meth_states, cutoff) {
  rle_states <- S4Vectors::Rle(meth_states)
  for (i in seq(2, S4Vectors::nrun(rle_states) - 1)){
    is_na <- is.na(S4Vectors::runValue(rle_states)[i])
    is_shorter <- S4Vectors::runLength(rle_states)[i] < cutoff
    has_same_flanking <- S4Vectors::runValue(rle_states)[i-1] == S4Vectors::runValue(rle_states)[i+1]
    if (is_na && (has_same_flanking || is_shorter)){
      meth_states[S4Vectors::start(rle_states)[i]:S4Vectors::end(rle_states)[i]] <-
        S4Vectors::runValue(rle_states)[i - 1]
    }
  }
  i <- 1
  is_na <- is.na(S4Vectors::runValue(rle_states)[i])
  is_shorter <- S4Vectors::runLength(rle_states)[i] < cutoff
  if (is_na && is_shorter){
    meth_states[S4Vectors::start(rle_states)[i]:S4Vectors::end(rle_states)[i]] <-
      S4Vectors::runValue(rle_states)[i + 1]
  }
  i <- S4Vectors::nrun(rle_states)
  is_na <- is.na(S4Vectors::runValue(rle_states)[i])
  is_shorter <- S4Vectors::runLength(rle_states)[i] < cutoff
  if (is_na && is_shorter){
    meth_states[S4Vectors::start(rle_states)[i]:S4Vectors::end(rle_states)[i]] <-
      S4Vectors::runValue(rle_states)[i - 1]
  }
  return(meth_states)
}
