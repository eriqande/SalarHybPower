
#' split data set into training and test and rank markers based on training set
#' 
#' 
#' We are going to have a single function that splits the same into training and test 
#' and also ranks the loci. 
#' @param dat a data frame that has the columns: id, variant, gene_copy, allele,
#' chrom, coord, pop, group.  
#' @param train_fract fraction of individuals from each group (and pop, too) to 
#' assign to the training set.
#' @param min_dist the number of base pairs below which you will not select a marker if
#' another one within min_dist on the same chromosome has already been selected.
#' @return a list with components `split_dat` and `ranked_markers`.
#' @export
split_and_rank <- function(dat, train_fract = 0.5, min_dist = 6e04) {
  A <- assign_train_test(dat, train_fract = train_fract) 
  R <- compute_epp(A)
  M <- marker_selection(R, min_dist = min_dist)
  list(split_dat = A, ranked_markers = M)
}
