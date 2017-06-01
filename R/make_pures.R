#' make pure individuals by segregating chunks of genome amongst founders
#'
#' Since founders are scrambled, this is not super meaningful, but it is still
#' a good way to simulate new genotypes.
#' internal function
#' @param Df a data frame of individuals with haplotypes in hap1 and hap2. If these are
#' founders then the haplotypes should have been scrambled first, using
#' \code{\link{scramble_founder_haplotypes}}.
#' @keywords internal
make_pures <- function(Df) {
  # basically, we split Df in half and then make F1s of them
  ids <- unique(Df$id)
  n <- length(ids)
  n2 <- floor(n/2)
  
  idsP <- sample(x = ids, size = n2 * 2, replace = FALSE)
  
  A <- Df %>%
    dplyr::filter(id %in% idsP[1:n2])
  B <- Df %>%
    dplyr::filter(id %in% idsP[(n2 + 1):(n2 * 2)])
  
  make_f1s(A, B)
}

