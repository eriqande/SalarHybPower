#' Compute the expected posterior probability of assignment of gene copies.
#'
#' See the notebooks for an explanation.
#' @param D a data frame like that returned by assign_train_test()
#' @return a data frame that includes variant, chrom, coord, allele then has
#' two columns of allele counts (farmed and wild) taken only from the "train"
#' individuals and then two columns of allele freqs (with the 1/2 correction
#' in there) and  epp computed from those allele freqs.
#' @keywords internal
compute_epp <- function(D) {
  D %>%
    filter(test_or_train == "train") %>%
    group_by(chrom, coord, variant, group, allele) %>%
    summarise(counts = n()) %>%
    tidyr::spread(data = ., key = group, value = counts, fill = 0) %>%
    rename(farmed_cnts = farmed,
           wild_cnts = wild) %>%
    group_by(chrom, coord, variant) %>%
    mutate(farmed_freq = (farmed_cnts + 1/2) / sum(farmed_cnts + 1/2),
           wild_freq = (wild_cnts + 1/2) / sum(wild_cnts + 1/2),
           epp_f = farmed_freq ^ 2 / (farmed_freq + wild_freq),  # epp_f and epp_w are just intermediate allele-specific
           epp_w = wild_freq ^ 2 / (farmed_freq + wild_freq),    # terms that make the calculation easier, here
           epp = (sum(epp_f) + sum(epp_w)) / 2)   %>%
    select(-epp_f, -epp_w)
}
