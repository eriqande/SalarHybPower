#' @param Wf a data frame of founders with scrambled haplotypes in hap1 and hap2
#' @param Ff a data frame of founders (from the other species/population) 
#' with scrambled haplotypes in hap1 and hap2
make_f1s <- function(Wf, Ff) {
  A <- Wf %>%
    count(id) %>%
    select(-n)
  B <- Ff %>%
    count(id) %>%
    select(-n)
  
  ML <- arrange_matings(A, B) # get the list of matings
  
  Wft <- Wf %>%
    select(id, variant, chrom, coord, hap1, hap2) %>%
    rename(id_A = id, `1` = hap1, `2` = hap2) %>%
    tidyr::gather(data = ., key = "gam_A", value = "hap1", `1`, `2`) %>%
    mutate(gam_A = as.numeric(gam_A))
  
  Fft <- Ff %>%
    select(id, variant, chrom, coord, hap1, hap2) %>%
    rename(id_B = id, `1` = hap1, `2` = hap2) %>%
    tidyr::gather(data = ., key = "gam_B", value = "hap2", `1`, `2`) %>%
    mutate(gam_B = as.numeric(gam_B))
  
  left_join(ML, Wft) %>%
    left_join(., Fft)
  
}

