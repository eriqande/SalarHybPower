#' Segregate the markers on two haploypes on each individual into two gametes.
#'
#' Used internally.
#' @param H a data frame like the one that comes out of scramble_founder_haplotypes().  It has
#' to have hap1 and hap2 columns.
#' @param Mb_recomb_prob the probability of recombination over 10^6 bp (one Mb).
#' By default this is set to 0.01 (i.e. a centi-Morgan).
#' @keywords internal
segregate_gametes <- function(H, Mb_recomb_prob = 0.01) {
  H %>%
    dplyr::arrange(id, chrom, coord) %>%
    dplyr::group_by(id, chrom) %>%
    dplyr::mutate(recom_prob = (coord - dplyr::lag(coord, 1)) * 0.01 / 10^6,
           recom_prob = ifelse(is.na(recom_prob), 0, recom_prob),
           start_value = sample(c(0,1), size = 1),
           cumul_xovers = start_value + cumsum(runif(n = n()) < recom_prob)) %>%
    # now, if cumum_xovers is even, gam1 is from hap1 and gam2 is from hap2.  Otherwise
    # gam1 is from hap2 and gam2 is from hap1
    dplyr::mutate(gam1 = ifelse(cumul_xovers %% 2 == 0, hap1, hap2),
           gam2 = ifelse(cumul_xovers %% 2 == 0, hap2, hap1)) %>%
    dplyr::select(-recom_prob, -start_value, -cumul_xovers) %>%
    dplyr::ungroup()
}
