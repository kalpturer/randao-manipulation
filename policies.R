# Example policies from the paper

tailmax <- function(len) {
  states <- expand.grid(x = -len:len, y = 0:len)
  colnames(states) <- c("val", "tail")
  states <- dplyr::filter(states, (tail == len & val <= 0) | (tail < len & val <= len - tail - 1))
  # sort states here
  prev_states <- states
  states <- dplyr::arrange(states, tail, val)

  states$prev_val <- c(NA, head(states[["val"]], -1))
  states$prev_tail <- c(NA, head(states[["tail"]], -1))
  return(states)
}

valuemax <- function(len) {
  states <- expand.grid(x = -len:len, y = 0:len)
  colnames(states) <- c("val", "tail")
  states <- dplyr::filter(states, (tail == len & val <= 0) | (tail < len & val <= len - tail - 1))

  score <- function(v, t) {
    return(v + t)
  }
  states$score <- mapply(score, states$val, states$tail)

  # sort states here
  prev_states <- states
  states <- dplyr::arrange(states, score, tail, val)

  states$prev_val <- c(NA, head(states[["val"]], -1))
  states$prev_tail <- c(NA, head(states[["tail"]], -1))
  return(states)
}