library(tidyverse)
library(ggridges)
library(superheat)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(modelr)

source("utilities.R")

results_dir <- "results_cauchy"

results_files <- dir(results_dir)
## results_files  <- results_files[grepl("2019-11-30", results_files)]
                                        #results_files  <- results_files[length(results_files)]

cols <- brewer.pal(8, "Set1")
cols_vec <- cols[c(1, 1, 1, 8, 8, 2, 2, 2, 3, 4, 5)]
lty_vec <- c("solid", "dashed", "dotted")[c(2, 3, 1, 1, 2, 2, 3, 1, 1, 1, 1)]

for(i in 1:length(results_files)) {

  file_name <- results_files[i]
  file_name
  file_path <- paste0(results_dir, "/", file_name)


  print(file_name)

  load(file_path)
  params <- stringr::str_match(file_name,
                               "results_n(\\d+)_p(\\d+)_coef(\\d*[\\.\\d]+)_escale(-?\\d+\\.?\\d*)_mscale(-?\\d+\\.?\\d*)_yalpha(\\d+)_talpha(\\d+)_estpropensity(TRUE|FALSE)_([ATE]+)_(-?\\d+\\.?\\d*)")
  n <- as.numeric(params[2])
  p <- as.numeric(params[3])
  coef <- as.numeric(params[4])
  escale <- as.numeric(params[5])
  mscale <- as.numeric(params[6])
  yalpha <- as.numeric(params[7])
  talpha <- as.numeric(params[8])
  estimated_propensity <- as.logical(params[9])
  estimand <- params[10]
  ab_dot_prod_true <- as.numeric(params[11])

  bias_var_plot  <- make_bias_var_plot(results_array, true_ate)

  plot_name <- sprintf("results_n%i_p%i_coef%.2f_escale%.2f_mscale%.2f_yalpha%i_talpha%i_estpropensity=%s_%s",
                       n, p, coef, escale, mscale, yalpha, talpha, estimated_propensity, estimand)

  ## rmse_plot + bias_plot + sd_plot +
  ##     plot_annotation(title = plot_title,
  ##                     subtitle = subtitle)
  if(estimated_propensity)
    ggsave(filename=sprintf("figs/%s.pdf", plot_name), bias_var_plot, width=7, height=3)
  else
    ggsave(filename=sprintf("figs_known_pscore/%s.pdf", plot_name), bias_var_plot, width=7, height=3)
}

