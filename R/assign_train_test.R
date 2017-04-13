#' assign to training and test
#'
#' Simple version that does a fixed fraction from each pop
#' @param dat data frame with columns id, variant, gene_copy, allele, 
#' chrom, coord, pop, and group.  
#' @param train_fract what fraction of each pop should be training (the rest is test)
#' @return A data frame just like dat, but with a column called "test_or_train" that
#' has entries of "test" or "train" as appropriate for each fish id.
assign_train_test <- function(dat, train_fract = 0.5) {
  d2 <- dat %>%
    distinct(id, pop, group) %>%
    group_by(group, pop) %>%
  mutate(test_or_train = sample(c(rep("train", ceiling(train_fract * n())),
                                  rep("test", n() - ceiling(train_fract * n()))
                                  )))
         
  left_join(dat, d2) 
}
