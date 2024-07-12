# The distributions of the count and the tail in addition to
# a helper function to precompute the distributions.

tail_dist_ar <- function(alpha, len) {
  tail_pdf = array(rep(NaN, len+1), len+1)

  for (i in 0:len) {
    if (i < len) {
      tail_pdf[i + 1] = ((1 - alpha) * (alpha ^ i))
    } else {
      tail_pdf[i + 1] = (alpha ^ len)
    }
  }

  return(tail_pdf)
}

count_dist_ar <- function(alpha, tail, len) {
  count_pdf = array(rep(NaN, (len - tail)), len - tail)

  for (i in 0:(len - tail - 1)) {
    count_pdf[i + 1] = choose(len - tail - 1, i) * (alpha ^ i) * ((1 - alpha) ^ (len - tail - 1 - i))
  }

  return(count_pdf)
}

precompute_dist <- function(alpha, len) {
  # precompute distributions
  t_pdf = tail_dist_ar(alpha, len)
  c_df = list()
  for (t in 0:(len-2)) {
    c_df[[t+1]] = count_dist_ar(alpha, t, len)
  }
  count_tail_dist_pdf <- function(v, t) {
    # P((C,T) = (v,t)) = P(T = t) * P(C = v | T = t)
    if (t == len || t == len-1) {
      if (v == 0) { return(t_pdf[t+1]) }
      else { return(0) }
    }
    if (v < 0 || v > len - t - 1 || t < 0 || t > len) {
      return(0)
    }

    return(t_pdf[t+1] * c_df[[t+1]][v+1])
  }
  return(count_tail_dist_pdf)
}