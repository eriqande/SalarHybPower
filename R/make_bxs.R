#' Make backcrosses by segregating chromosome segments
#'
#' blab
#' @inheritParams make_f1s
#' @param idstr string that gets prepended to the number of each individual
#' @keywords internal
make_bxs <- function(Wf, Ff, idstr = "BX") {
  A <- Wf %>%
    count(id) %>%
    select(-n)
  B <- Ff %>%
    count(id) %>%
    select(-n)

  M <- min(floor(nrow(A)/3), nrow(B))

  Aids <- sample(A$id) # permute these
  Bids <- sample(B$id) # permute these too

  Ai1 <- Aids[1:M]  # ids of A's for F1 creation
  Bi1 <- Bids[1:M]  # ids of B's for F1 creation

  Ai2 <- Aids[(M + 1):(3 * M)]  # ids of A's for making backcrosses

  # make the F1s and segregate em
  F1 <- make_f1s(Wf %>% filter(id %in% Ai1),
                 Ff %>% filter(id %in% Bi1)) %>%
    segregate_gametes()

  # now, make F1's between the F1s and the Ai2's
  BX <- make_f1s(Wf %>% filter(id %in% Ai2), F1) %>%
    mutate(id = stringr::str_replace(id, "^F1", idstr))

}
