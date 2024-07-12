# An implementation of the policy iteration algorithm as defined in the paper

stationary_distribution <- function(transition, len) {
  # https://stephens999.github.io/fiveMinuteStats/markov_chains_discrete_stationary_dist.html
  A <- t(diag(rep(1,len+1))-transition)
  b <- rep(0,len+1)
  AA <- rbind(A,rep(1,len+1))
  bb <- c(b,1)
  pi_lineq <- t(solve(t(AA)%*%AA,t(AA)%*%bb))
  return(pi_lineq)
}

policy_eval <- function(alpha, len, states, count_tail_dist_pdf) {
  # transition_probabilities
  P = matrix(data = NA,
             nrow = len + 1,
             ncol = len + 1)
  # exprected reward from state
  R_e = rep(NaN, len+1)
  # values for states
  values = rep(NaN, len+1)

  maxv = array(rep(NaN, (len+1)*(len+1)*(2*len+1)), c(len+1, len+1, (2*len+1)))
  for (t_src in (0:len)) {
    # precompute cumulative sum
    debug_println(paste("Computing maxpair cdf:", t_src, "out of", len))
    csum = array(rep(0, (t_src+1)), c(t_src+1))
    for(r in 1:nrow(states)) {
      v = states[r, "val"]
      t_dst = states[r, "tail"]
      prod = 1
      for (i in 0:t_src) {
        csum[i+1] = csum[i+1] + count_tail_dist_pdf(v + i, t_dst)
        if (csum[i+1] <= 1) {
             prod = prod * (csum[i+1] ^ (choose(t_src, i)))
        }
      }
      maxv[t_src + 1, t_dst + 1, v + (len + 1)] = prod
    }
  }
  debug_print("maxv_last: ")
  debug_println(maxv[,states[nrow(states), "tail"]+1,states[nrow(states), "val"]+(len+1)])
  debug_print("maxv_last error: ")
  debug_println(max(abs(1 - maxv[,states[nrow(states), "tail"]+1,states[nrow(states), "val"]+(len+1)])))

  for (t_src in (0:len)) {
    debug_println(paste("Computing transition probabilities:", t_src, "out of", len))
    R_ans = 0
    for (t_dst in (0:len)) {
      ansP = 0
      vrange = ifelse(t_dst == len, 0, len-t_dst-1)
      for (v in (-t_src):vrange) {
        row = dplyr::filter(states,  tail %in% t_dst, val %in% v)
        if (is.na(row$prev_tail) && is.na(row$prev_val)) {
          ansP = ansP + maxv[t_src + 1, t_dst + 1, v + (len + 1)]
          R_ans = R_ans + (maxv[t_src + 1, t_dst + 1, v + (len + 1)] * (t_dst + v))
        } else {
          ansP = ansP + (maxv[t_src + 1, t_dst + 1, v + (len + 1)] - maxv[t_src + 1, row$prev_tail + 1, row$prev_val + (len + 1)])
          R_ans = R_ans + ((maxv[t_src + 1, t_dst + 1, v + (len + 1)] - maxv[t_src + 1, row$prev_tail + 1, row$prev_val + (len + 1)])*(t_dst + v))
        }
      }
      P[t_src+1, t_dst+1] = ansP
    }
    R_e[t_src+1] = R_ans
  }
  pi = stationary_distribution(P, len)
  reward <- sum(R_e * pi)

  debug_print("rowSums(P) = ")
  debug_println(rowSums(P))
  debug_print("rowSums(P) error: ")
  debug_println(max(abs(1 - rowSums(P))))

  A = diag(len+1) - P
  A = A[,2:(len+1)]
  b = R_e - reward
  values <- qr.solve(A, b, tol = 0)
  values <- c(0, values)

  return(list(reward, values))
}

policy_improvement <- function(alpha, len, values, reward) {
  # (V,T) pairs generate where V is the nontail value
  states <- expand.grid(x = -len:len, y = 0:len)
  colnames(states) <- c("val", "tail")
  # filter invalid rows
  states <- dplyr::filter(states, (tail == len & val <= 0) | (tail < len & val <= len - tail - 1))

  bellman <- function(v, t) {
    return(values[t+1] + v + t)
  }
  states$score <- mapply(bellman, states$val, states$tail)

  # sort states here
  prev_states <- states
  states <- dplyr::arrange(states, score, tail, val)

  states$prev_val <- c(NA, head(states[["val"]], -1))
  states$prev_tail <- c(NA, head(states[["tail"]], -1))
  return(states)
}


policy_iteration <- function(alpha, len, max_iter) {
  count_tail_dist_pdf <- precompute_dist(alpha, len)

  # initialize random policy
  states <- expand.grid(x = -len:len, y = 0:len)
  colnames(states) <- c("val", "tail")
  states <- dplyr::filter(states, (tail == len & val <= 0) | (tail < len & val <= len - tail - 1))
  states <- states[sample(1:nrow(states)), ]
  prev_states <- states
  states$prev_val <- c(NA, head(states[["val"]], -1))
  states$prev_tail <- c(NA, head(states[["tail"]], -1))

  # run policy iteration
  for (i in (1:max_iter)) {
    debug_println(i)
    l = policy_eval(alpha, len, states, count_tail_dist_pdf)
    reward = l[[1]]
    values = l[[2]]

    new_states = policy_improvement(alpha, len, values, reward)
    debug_println(reward)
    prev_states = states
    states = new_states

    if (identical(states, prev_states)) {
      debug_println("CONVERGED")
      debug_print("final reward: ")
      debug_println(reward)
      debug_print("final values: ")
      debug_println(values)
      return(list(reward, values, new_states))
    }
  }
}