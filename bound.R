# Script to compute the bounds for the appendix of the full paper
library(ggplot2)
set.seed(1234)
options(digits = 17)

# parameters
LEN = 32
alphas = seq(0.01, 0.49, 0.01)
debug = TRUE

# helper functions
debug_print <- function(str) {
  if (debug) { cat(str) }
}
debug_println <- function(str) {
  if (debug) { print(str) }
}
map_to_index <- function(x, y) {
  return(x * (LEN + 1) + y + 1)
}

source('./policies.R')
source('./distributions.R')
source('./policy_iteration.R')


bounds = c()
for (alpha in alphas) {
  print(paste("alpha =", alpha))

  # compute transition matrix
  count_tail_dist_pdf <- precompute_dist(alpha, LEN)
  tm_results = policy_eval(alpha, LEN, tailmax(LEN), count_tail_dist_pdf)
  P = tm_results[[3]]


  ## Encode system of linear equations from the appendix
  # helper functions

  # E[H(t,L)] \leq bound
  tree_bound <- function(t) {
    return(1 + (1 / (1 - ((2 * alpha)**LEN) * (1 + ((2 * alpha)**LEN)*(1 - alpha**LEN)))))
  }

  A = matrix(0L,
             nrow = (LEN + 1)*(LEN + 1),
             ncol = (LEN + 1)*(LEN + 1))
  b = rep(0, (LEN + 1)*(LEN + 1))

  i = 1
  for (t in 0:LEN) {
    A[i, map_to_index(t,0)] = 1
    b[i] = 0
    i = i + 1
  }
  for (t in 0:LEN) {
    A[i, map_to_index(t,LEN)] = 1
    A[i, map_to_index(LEN-1,LEN-1)] = -1
    b[i] =  tree_bound(t)
    i = i + 1
  }
  for (t in 0:LEN) {
    for (tt in 1:(LEN-1)) {
      A[i, map_to_index(t,tt)] = 1
      for (ttt in 1:LEN)  {
        A[i, map_to_index(tt,ttt)] = -(P[t+1, ttt+1])
        if (t == tt && tt == ttt) {
          A[i, map_to_index(tt,ttt)] = A[i, map_to_index(tt,ttt)] + 1
        }
      }
      b[i] = 1
      i = i + 1
    }
  }

  maxval <- max(qr.solve(A, b))
  print(paste("max_bound: ", maxval))
  bounds <- c(bounds, maxval)
}

df = data.frame('alpha' = alphas,
                'bound' = bounds)

# save results
data.table::fwrite(file = "bounds.csv", df, row.names = FALSE, quote = FALSE, sep = ",")