#' assign individuals to training and test sets.
#'
#' Simple version that does a fixed fraction from each pop
#' @param dat data frame with columns id, variant, gene_copy, allele,
#' chrom, coord, pop, and group.
#' @param train_fract what fraction of each pop should be training (the rest is test)
#' @return A data frame just like dat, but with a column called "test_or_train" that
#' has entries of "test" or "train" as appropriate for each fish id.
#' @keywords internal
assign_train_test <- function(dat, train_fract = 0.5) {
  d2 <- dat %>%
    dplyr::distinct(id, pop, group) %>%
    dplyr::group_by(group, pop) %>%
  dplyr::mutate(test_or_train = sample(c(rep("train", ceiling(train_fract * n())),
                                  rep("test", n() - ceiling(train_fract * n()))
                                  )))

  dplyr::left_join(dat, d2)
}
