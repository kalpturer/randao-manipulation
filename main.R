library(dplyr)
library(ggplot2)
library(data.table)
library(parallel)
library(parabar)
set.seed(1234)
options(digits = 17)

# parameters
debug = FALSE
LEN = 8
alpha_from = 0
alpha_to   = 1
alpha_by   = 0.25
run_in_parallel = FALSE

# helper functions
debug_print <- function(str) {
  if (debug) { cat(str) }
}
debug_println <- function(str) {
  if (debug) { print(str) }
}

# import functions
source('./policies.R')
source('./distributions.R')
source('./policy_iteration.R')

# run policy iteration and generate results
alphas = seq(alpha_from, alpha_to, alpha_by)
parameters_df = data.frame('alpha_from' = alpha_from,
                           'alpha_to' = alpha_to,
                           'alpha_by' = alpha_by,
                           'LEN' = LEN)
folder = paste('./results_L', LEN, '/', sep = '')
if (!file.exists(folder)){
  dir.create(folder)
}
data.table::fwrite(file = paste(folder, 'parameters.csv', sep=''), parameters_df, quote = FALSE, sep = ",")


plot_results = data.frame('alpha' = alphas,
                          'optimal' = NA,
                          'honest' = alphas * LEN,
                          'tailmax' = NA,
                          'valuemax' = NA)


run <- function(i) {
  print(paste("alpha =", alphas[i]))
  opt_reward = NA
  tm = list(NA)
  vm = list(NA)
  m_error = -1
  p_error = -1

  try ({
    result_list = policy_iteration(alphas[i], LEN, 100)
    opt_reward = result_list[[1]]
    values = result_list[[2]]
    policy = result_list[[3]]
    policy = policy[,c('tail', 'val' , 'score')]
  })

  # precompute distributions
  count_tail_dist_pdf <- precompute_dist(alphas[i], LEN)

  try ({
    tm_results = policy_eval(alphas[i], LEN, tailmax(LEN), count_tail_dist_pdf)
    tm = tm_results[[1]]
  })
  try ({
    vm_results = policy_eval(alphas[i], LEN, valuemax(LEN), count_tail_dist_pdf)
    vm = vm_results[[1]]
  })

  return(
    list(
      alphas[i],
      opt_reward,
      tm[[1]],
      vm[[1]],
      alphas[i] * LEN
    )
  )
}

# run analysis
if (run_in_parallel) {
  set_option("progress_track", TRUE)
  be <- start_backend(cores = detectCores() - 1, cluster_type = "fork", backend_type = "async")
  export(be, ls(environment()), environment())
  results_list <- par_sapply(backend = be, 1:length(alphas), run)
  plot_results <- as.data.frame(t(matrix(unlist(results_list), ncol = length(alphas))))
  colnames(plot_results) <- c("alpha", "optimal", "tailmax", "valuemax", "honest")
  stop_backend(be)

  plot_results$tailmax = as.numeric(plot_results$tailmax)
  plot_results$valuemax = as.numeric(plot_results$valuemax)
} else {
  for (i in 1:length(alphas)) {
    res <- run(i)
    plot_results$optimal[i] = res[[2]]
    plot_results$tailmax[i] = res[[3]]
    plot_results$valuemax[i] = res[[4]]
    plot_results$honest[i] = res[[5]]
  }
}

# save results
data.table::fwrite(file = paste(folder, "results.csv", sep=""), plot_results, row.names = FALSE, quote = FALSE, sep = ",")

# plot results
plot_results$tailmax = as.numeric(plot_results$tailmax)
plot_results$valuemax = as.numeric(plot_results$valuemax)
cs <- c("honest" = "blue",
        "optimal" = "red",
        "tailmax" = "green",
        "valuemax" = "purple",
        "alpha" = "black")
g <- ggplot(data = plot_results, aes(x = alpha)) +
  geom_line(data = plot_results, aes(x = alpha, y = honest/LEN,   color = "honest")) +
  geom_line(data = plot_results, aes(x = alpha, y = tailmax/LEN,  color = "tailmax")) +
  geom_line(data = plot_results, aes(x = alpha, y = valuemax/LEN, color = "valuemax")) +
  geom_line(data = plot_results, aes(x = alpha, y = optimal/LEN,  color = "optimal")) +
  labs(x = "alpha",
       y = "reward",
       color = "Legend") +
  scale_color_manual(values = cs)
ggsave(
  paste(folder, "results_", LEN, ".pdf", sep=""),
  plot = g,
  device = "pdf",
  scale = 1,
  dpi = 600,
)


