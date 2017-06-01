#' write a new_hybs file
#'
#' Takes care of the preamble and everthing.
#' @param G a data frame like that prepared by nh_tablify, or perhaps several of those
#' that have been bind_rows()-ed together.
#' @param path  the file path you want to write to.
#' @keywords internal
write_nh <- function(G, path = "nh_data.txt") {
  # first write the preamble
  cat("NumIndivs ", nrow(G), "\nNumLoci ", ncol(G) - 2, "\nDigits 1\nFormat NonLumped\n\n", file = path)

  # then the locus names
  cat(c("LocusNames  ", names(G)[-(1:2)]), sep = " ", eol = "\n", file = path, append = TRUE)

  # then format the genotype frame and write it
  G %>%
    mutate(id = paste("n", id, sep = " ")) %>%
    ungroup() %>%
    mutate(idx = 1:n()) %>%
    select(idx, everything()) %>%
    write.table(file = path, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "   ", eol = "\n")

}
