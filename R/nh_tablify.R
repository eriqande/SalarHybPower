#' turn the genotypes of individuals into rows for a new hybrids file
#'
#' This also makes both gene copies missing if only one of them was.  This happens in the
#' simulated hybrids, so obviously must be dealt with.
#' @param D a tidy data frame with at least the columns: id, variant, gene_copy, allele
#' @param V a vector of variant names giving which variants to use and the order they are
#' to be in.
#' @param opt_str  The options string to be added to each individual (like "z0s")
#' @param id_prepend this string will be prepended to the fish id (i.e. use it for "train_", or "test_")
#' @keywords internal
nh_tablify <- function(D, V, opt_str = "", id_prepend = "") {
  D1 <- D %>%
    filter(variant %in% V)

  common_loc <- V[V %in% D1$variant]  # these are the loci from V that are in D, in V-order

  # now spread the gene copies and then the variants
  D1 %>%
    select(id, variant, gene_copy, allele) %>%
    tidyr::spread(data = ., key = gene_copy, value = allele) %>%
    mutate(gc1 = as.integer(ifelse(is.na(`1`) | is.na(`2`), -1, `1`)),  # make both -1 if either of them is -1
           gc2 = as.integer(ifelse(is.na(`2`) | is.na(`1`), -1, `2`))) %>%
    select(-`1`, -`2`) %>%
    mutate(geno = paste(gc1, gc2, sep = " ")) %>%
    select(-gc1, -gc2) %>%
    tidyr::spread(data = ., key = variant, value = geno, fill = "-1 -1") %>%
    .[, c("id", common_loc)] %>%
    mutate(option = opt_str) %>%
    mutate(id = paste(id_prepend, id, sep = "")) %>%
    select(id, option, everything())

}
