#' return list of  who mates with whom to create hybrid offspring.  This will create F1's between 
#' the A and the B group.
#' @param A a data frame with id, pop, and group.  Should be exclusively from one of the groups (i.e. wild)
#' @param B a data frame with id and pop, and group.  Should be from the other group (i.e. farmed)
#' @return a data frame with columns id, popA, popB, idA gamA idB gamB that tells us which gametes from
#' which individuals each new individual is to be made of.
arrange_matings <- function(A, B) {
  
  n <-  min(nrow(A), nrow(B)) 
  Atib <- tibble(id = rep(sample(A$id, size = n, replace = FALSE), each = 2))  %>%
    mutate(gam = rep(c(1,2), length.out = n())) %>%
    sample_n(size = nrow(.), replace = FALSE)
  Btib <- tibble(id = rep(sample(B$id, size = n, replace = FALSE), each = 2))  %>%
    mutate(gam = rep(c(1,2), length.out = n())) %>%
    sample_n(size = nrow(.), replace = FALSE)
  
  
  Atib2 <- left_join(Atib, A)
  names(Atib2) <- paste(names(Atib2), "A", sep = "_")
  Btib2 <- left_join(Btib, B)
  names(Btib2) <- paste(names(Btib2), "B", sep = "_")
  
  F1s <- bind_cols(Atib2, Btib2) %>%
    mutate(id = paste("F1", 1:n(), sep = "_")) %>%
    select(id, everything())
  
  return(F1s)
}

