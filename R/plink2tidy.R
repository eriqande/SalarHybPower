#' Convert a plink file to a tidy data frame
#'
#' This  might be useful to users.
#' @param ped path to a space delimited plink ped file
#' @param map path to a tab delimited plink map file
#' @export
plink2tidy <- function(ped, map) {
  m <- read_delim(map, delim = "\t", col_names = FALSE, progress = FALSE) %>%
    setNames(c("chrom", "variant", "dummy", "coord"))

  p <- read_delim(ped, delim = " ", col_names = FALSE, progress = FALSE, na = "0") %>%
    setNames(c("id", "fid", "d1", "d2", "d3", "d4", paste(rep(m$variant, each = 2), c(1,2), sep = " ")))

  # now we gather p and drop the columns that are not relevant here
  # and then split the gene copies up
  pt <- p %>%
    select(-(fid:d4)) %>%
    tidyr::gather(data = ., key = "snp_and_gc", value = "allele", -id) %>%
    tidyr::separate(data = ., col = snp_and_gc, into = c("variant", "gene_copy"), sep = " ", convert = TRUE)


  # now bind the marker information on there and make chrom a factor so it sorts in the right order
  mp <- left_join(pt, m %>% select(-dummy)) %>%
    mutate(chrom = factor(chrom, levels = unique(m$chrom))) %>%
    arrange(chrom, coord, id, gene_copy)

  # return that
  mp

}
