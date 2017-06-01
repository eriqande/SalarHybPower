#' Select most informative markers on chromosomes
#'
#' Does not select markers that are within min_dist base pairs
#' of a previously selected marker.
#' @param R a data frame like that returned by compute_epp()
#' @return returns a data frame with chrom, coord, variant, epp,
#' and a new column "selectable" which is TRUE if it is one that should
#' be selected, and FALSE if it is too close to another, previously selected
#' marker.  It comes back sorted, descending on epp.
#' @keywords internal
marker_selection <- function(R, min_dist = 6e04) {
  # first we need to whittle R down to only 1 row per variant, and just keep the epp
  # and we want to arrange on chrom then epp
  R1 <- R %>%
    group_by(chrom, coord, variant) %>%
    summarise(epp = first(epp)) %>%
    ungroup() %>%
    arrange(chrom, desc(epp))

  R1 %>%
    arrange(chrom, desc(epp)) %>%
    group_by(chrom) %>%
    mutate(selectable = ranked_selector(coord, min_dist)) %>%
    ungroup() %>%
    arrange(desc(epp))
}



#' an internal function to be applied to each chromosome
#'
#' this is a very internal function.  Users should never see this.
#' @param co a vector of base pair coordinates ranked descending on epp
#' @param md the min_dist allowed between a new SNP and any other selected ones
#' @keywords internal
ranked_selector <- function(bp, md) {
  N <- length(bp)
  ret <- rep(NA, N)
  ret[1] <- TRUE  # initialize by always choosing the first SNP
  n <- 1  # the number chosen
  B <- bp[1]

  if (N == 1) return(ret)

  for (i in 2:N) {
    if (all(abs(bp[i] - B) > md)) {
      ret[i] <- TRUE
      n <- n + 1
      B[n] <- bp[i]

    } else {
      ret[i] <- FALSE
    }
  }
  ret
}

