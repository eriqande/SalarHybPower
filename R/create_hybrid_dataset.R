
#' Create a simulated data set for newhybrids using dplyr::selected/ranked markers
#' 
#' We do this so that all the training individuals will end up with `z0s` or `z1s` in the data
#' set, and then we make as many of the others as we can.  Since we are going to be just computing 
#' the scaled likelihoods for all these guys, it is fine to just create everyone to be of one type
#' (i.e. Pure, F1, F2, Bx, etc.)  This let's you make sure that backcrosses are all within the same
#' population.
#' @param SAR the output of split_and_rank()
#' @param wild_pop a vector of the names of the populations of wild fish to simulate from. It treats them
#' as all being from the same population if they are included together.  This is not a 
#' problem for F1's, but it is likely unrealistic for F2s and backcrosses. Unless the 
#' F1 hybrids have no site fidelity (i.e. they stray freely...)
#' @param hyb_cat the hybrid category desired 
#' @param L the number of loci to choose
#' @param dir the path of the directory to create to put the result into.  It will 
#' write the result into a file called \code{file_name}.
#' @param file_name the name of the output file.  Defaults to "nh_data.txt".
#' @export
create_hybrid_dataset <- function(SAR, wild_pop, hyb_cat, L, dir = "hyb_dir", file_name = "nh_data.txt") {
  
  # Get the L markers. store their names in a vector called V (for variants)
  V <- SAR$ranked_markers %>%
    dplyr::filter(selectable == TRUE & cumsum(selectable) <= L) %>%
    .$variant  
  
  # break the data into training and test individuals
  Train <- SAR$split_dat %>%
    dplyr::filter(test_or_train == "train")
  Test <- SAR$split_dat %>%
    dplyr::filter(test_or_train == "test")
  
  # make the rows that the Train indivs will be in the newhybrids data set.
  # note that "farmed" will always be species 0. 
  trainNH <- list(
    Train %>%
      dplyr::filter(group == "farmed") %>%
      nh_tablify(., V, opt_str = "z0s", id_prepend = "train_"),
    Train %>%
      dplyr::filter(group == "wild") %>%
      nh_tablify(., V, opt_str = "z1s", id_prepend = "train_")
  ) %>%
    dplyr::bind_rows()
  
  # now make the data frames that represent the wild and farmed founders and 
  # scramble their haplotypes up
  Wf <- Test %>%
    dplyr::filter(group == "wild") %>%
    dplyr::filter(pop %in% wild_pop) %>%
    scramble_founder_haplotypes(., V)
  
  Ff <- Test %>%
    dplyr::filter(group == "farmed") %>%
    scramble_founder_haplotypes(., V)
  
  
  # now, simulate whatever type of hybrid was requested
  Hybs <- switch(hyb_cat,
                 PureW = make_pures(Wf),
                 PureF = make_pures(Ff), 
                 F1 = make_f1s(Wf, Ff),
                 F2 = make_f2s(Wf, Ff),
                 BX = make_bxs(Wf, Ff),
                 stop("unknown hyb_cat: ", hyb_cat)
  )
  
  # and prepare each of those simulated individuals as a row in new-hybs file
  hybsNH <- Hybs %>%
    dplyr::rename(`1` = hap1, `2` = hap2) %>%
    tidyr::gather(data = ., key = "gene_copy", value = "allele", `1`, `2`) %>%
    dplyr::mutate(gene_copy = as.integer(gene_copy)) %>%
    nh_tablify(., V, opt_str = "", id_prepend = "")
  
  # now, make the directory to write the result into
  dir.create(dir)
  dplyr::bind_rows(trainNH, hybsNH) %>%
    write_nh(., path = file.path(dir, file_name))
  
}
