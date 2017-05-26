make_f2s <- function(Wf, Ff) {
  F1 <- make_f1s(Wf, Ff) %>%
    segregate_gametes()
  
  # now, make two new data frames that we will combine using make_F1s into F2s
  ids <- sample(unique(F1$id))  # permute all the ids around
  if (length(ids) %% 2 == 1) {  # if there is an odd number of F1s, then drop the first one
    ids <- ids[-1]
  }
  idA <- ids[1:(length(ids) / 2)]
  idB <- ids[(1 + length(ids) / 2):length(ids)]
  
  # now make new Wf and new Ff and run them back through make_F1s
  F1_clean <- F1 %>%
    select(-ends_with("_A"), -ends_with("_B"))
  
  newWf <- F1_clean %>% filter(id %in% idA)
  newFf <- F1_clean %>% filter(id %in% idB)
  
  F2 <- make_f1s(newWf, newFf) %>%
    mutate(id = stringr::str_replace(id, "^F1", "F2"))
  
}
