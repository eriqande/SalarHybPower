#' Randomly assign (no-LD) markers to one chromosome or another in  founder individuals
#'
#' blah blah
#' @param D a tidy data frame with at least the columns: id, variant, gene_copy, allele, chrom, coord, pop
#' @param V a vector of variant names giving which variants to use and the order they are
#' to be in.
#' @return A data frame with all the identifiers as before, but now the alleles are in two
#' columns: hap1 and hap2, which are the alleles as if they occurred together on different
#' homologous chromosomes.  This gets sorted by id, chrom, and coord.
#' @keywords internal
scramble_founder_haplotypes <- function(D, V) {
  # make one column for each gene copy
  D1 <- D %>%
    dplyr::ungroup() %>%
    dplyr::filter(variant %in% V) %>%
    tidyr::spread(data = ., key = gene_copy, value = allele) %>%
    dplyr::mutate(segind = sample(x = c(T,F), size = n(), replace = TRUE),
           hap1 = ifelse(segind, `1`, `2`),
           hap2 = ifelse(segind, `2`, `1`)) %>%
    dplyr::select(-`1`, -`2`, -segind) %>%
    dplyr::arrange(id, chrom, coord)

}
