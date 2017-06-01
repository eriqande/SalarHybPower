#' Mate two collections of individuals together to make as many offspring as possible
#'
#' This does sampling without replacement from the genes of the individuals.
#' Basically, it takes individuals that have haplotypes in hap1 and hap2, and it
#' segregates two gametes for each individual (with recombination) and then those
#' get combined with gametes from other individuals.  If there are W indivs in
#' Wf and F indivs in Ff, this function will return 2 * min(W, F) offspring from mating between
#' min(W, F) randomly selected individuals in Wf and Ff.
#' @param Wf a data frame of individuals with haplotypes in hap1 and hap2. If these are
#' founders then the haplotypes should have been scrambled first, using
#' \code{\link{scramble_founder_haplotypes}}.
#' @param Ff a data frame of individuals (from the other species/population)
#' with haplotypes in hap1 and hap2. If these are
#' founders then the haplotypes should have been scrambled first, using
#' \code{\link{scramble_founder_haplotypes}}.
#' @keywords internal
make_f1s <- function(Wf, Ff) {
  A <- Wf %>%
    count(id) %>%
    select(-n)
  B <- Ff %>%
    count(id) %>%
    select(-n)

  ML <- arrange_matings(A, B) # get the list of matings

  Wft <- Wf %>%
    select(id, variant, chrom, coord, hap1, hap2) %>%
    segregate_gametes() %>%
    rename(id_A = id, `1` = gam1, `2` = gam2) %>%
    select(-hap1, -hap2) %>%
    tidyr::gather(data = ., key = "gam_A", value = "hap1", `1`, `2`) %>%
    mutate(gam_A = as.numeric(gam_A))

  Fft <- Ff %>%
    select(id, variant, chrom, coord, hap1, hap2) %>%
    segregate_gametes() %>%
    rename(id_B = id, `1` = gam1, `2` = gam2) %>%
    select(-hap1, -hap2) %>%
    tidyr::gather(data = ., key = "gam_B", value = "hap2", `1`, `2`) %>%
    mutate(gam_B = as.numeric(gam_B))

  left_join(ML, Wft) %>%
    left_join(., Fft)

}

