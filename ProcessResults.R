library(tidyverse)
library(ggridges)
library(superheat)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(modelr)

results_dir <- "results"

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
    
    rmse_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))
    bias_mat <- abs(apply(results_array - true_ate, c(2, 3), function(x) mean(x, na.rm=TRUE)))
    var_mat <- rmse_mat^2 - bias_mat^2

    iters <- dim(results_array)[1]
    se_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) sd(x, na.rm=TRUE))/iters)
    se_tib  <- as_tibble(t(se_mat))
    se_tib$W <- as.numeric(colnames(rmse_mat))
    
    
    tib <- as_tibble(t(rmse_mat))

    tib$W <- as.numeric(colnames(rmse_mat))

    tib2  <- inner_join(tib %>% gather(key=Type, value=RMSE, -W),
                        se_tib  %>% gather(key=Type, value=se, -W),
                        by=c("W", "Type"))
    
    rmse_plot <- tib %>% gather(key=Type, value=RMSE, -W) %>%
        ggplot() + geom_line(aes(x=W, y=RMSE, col=Type, linetype=Type)) +
        theme_bw() + theme(legend.position="none") +
        scale_color_manual(values=cols_vec) +
        scale_linetype_manual(values=lty_vec) +
        xlab(expression(w[2]))
    
    rmse_plot <- tib2 %>%  ggplot() + geom_line(aes(x=W, y=RMSE, col=Type, linetype=Type)) +
        geom_linerange(aes(x=W, ymin=sqrt(RMSE^2-1.96*se), ymax=sqrt(RMSE^2+1.96*se), col=Type, linetype=Type), size=0.2) +
        theme_bw() + theme(legend.position="none") +
        scale_color_manual(values=cols_vec) +
        scale_linetype_manual(values=lty_vec) +
        xlab(expression(w[2]))

    tib <- as_tibble(t(bias_mat))
    tib$W <- as.numeric(colnames(bias_mat))
    bias_plot <- tib %>% gather(key=Type, value=Bias, -W) %>%
        ggplot() + geom_line(aes(x=W, y=Bias, col=Type, linetype=Type)) +
        theme_bw() + theme(legend.position="none") + geom_hline(yintercept=0, linetype=2)  +
        scale_color_manual(values=cols_vec) +
        scale_linetype_manual(values=lty_vec) +
      xlab(expression(w[2])) +
      ylab("Absolute Bias")

    tib <- as_tibble(t(sqrt(var_mat)))
    tib$W <- as.numeric(colnames(bias_mat))
    sd_plot <- tib %>% gather(key=Type, value=SD, -W) %>%
        ggplot() + geom_line(aes(x=W, y=SD, col=Type, linetype=Type)) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=cols_vec) +
        scale_linetype_manual(values=lty_vec) +
        xlab(expression(w[2]))

  if(!is.null(ab_dot_prod_true))
    ab_dot_prod <- ab_dot_prod_true

        

    if(yalpha == 0) 
        plot_title <- "Ridge Regression"
    if(yalpha == 1)
        plot_title <- "Lasso Regression"
    

    if(estimated_propensity)
        plot_title <- paste(plot_title, "Estimated Propensity Score", sep=", ")
    if(!estimated_propensity)
        plot_title <- paste(plot_title, "Known Propensity Score", sep=", ")

    ## Plot
    
    plot_name <- sprintf("results_n%i_p%i_coef%.2f_escale%.2f_mscale%.2f_yalpha%i_talpha%i_estpropensity=%s_%s",
                         n, p, coef, escale, mscale, yalpha, talpha, estimated_propensity, estimand)

    ab_dot_prod_round <- round(ab_dot_prod, 2)
    subtitle <- bquote("n=" ~.(n) ~ ", p=" ~.(p) ~ "," ~ s[T] == .(escale) ~ "," ~ s[Y] == .(mscale) ~ "," ~ alpha^T ~ beta ~ "=" ~ .(ab_dot_prod_round))

    rmse_plot + bias_plot + sd_plot +
        plot_annotation(title = subtitle)



    ## rmse_plot + bias_plot + sd_plot +
    ##     plot_annotation(title = plot_title,
    ##                     subtitle = subtitle)
    if(estimated_propensity)
      ggsave(filename=sprintf("figs_icml/%s.pdf", plot_name), width=7, height=3)
    else
      ggsave(filename=sprintf("figs_icml_known_pscore/%s.pdf", plot_name), width=7, height=3)
}


## results_dir <- "results"
results_files <- dir(results_dir)


## file_matches <- stringr::str_extract(results_files,
##                                      sprintf("results_n%i_p%i_coef0_escale(-?\\d+\\.?\\d*)_mscale(\\d+\\.?\\d*)_yalpha[0-1]+_estpropensity%s", n, p, estpropensity))
## file_matches <- results_files[!is.na(file_matches)]
## file_matches <- file_matches[-c(6, 8, 10, 12, 18, 20, 22, 24)]

file_matches <- results_files[-c(1:6)]
file_matches


rmse_table <- matrix(nrow=10, ncol=length(file_matches))
table_colnames <- c()
count <- 1

escale_vec <- c()
mscale_vec <- c()
yalpha_vec <- c()


for(fn in file_matches){

    file_path <- paste(results_dir, fn, sep="/")

    load(file_path)
    params <- stringr::str_match(fn, "results_n(\\d+)_p(\\d+)_coef([0-9])+_escale(-?\\d+\\.?\\d*)_mscale(-?\\d+\\.?\\d*)_yalpha(\\d+)_estpropensity(TRUE|FALSE)")

    n <- as.numeric(params[2])
    p <- as.numeric(params[3])
    coef <- as.numeric(params[4])
    escale <- as.numeric(params[5])
    mscale <- as.numeric(params[6])
    yalpha <- as.numeric(params[7])
    estimated_propensity <- as.logical(params[8])

    rmse_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))

    rmse_table[, count] <- c(rmse_mat[c("Naive", "Regression", "IPW", "AIPW"), 1],
                             rmse_mat[c("IPW_d_oracle", "IPW_d", "AIPW_d"), "0"],
                             rmse_mat[c("IPW_d_oracle", "IPW_d", "AIPW_d"), "-1"])
    
    table_colnames <-  c(table_colnames, paste(coef, escale, mscale, yalpha, estimated_propensity, sep=","))

    escale_vec <- c(escale_vec, escale)
    mscale_vec <- c(mscale_vec, mscale)
    yalpha_vec <- c(yalpha_vec, yalpha)
    
    count <- count + 1
}

rownames(rmse_table) <- c("Naive", "Regression", "IPW", "AIPW", "IPW-oracle-d(0)", "IPW-d(0)", "AIPW-d(0)", "IPW-oracle-d(1)", "IPW-d(-1)", "AIPW-d(-1)")

rmse_table <- apply(rmse_table, 2, function(x) {
    x <- round(x, 2)
    min_index <- which.min(x)
    x[min_index] <- paste0("textbf{", x[min_index], "}")
    x
})

## overlap_colname <- c("High Overlap"=8, "Low Overlap"=8)
## snr_colname <- rep(c("High SNR"=4, "Low SNR"=4), 2)
## yalpha_colname <- rep(2, 8)
## names(yalpha_colname) <- rep(c("R", "L"), 4)



dimension_colname <- c("Lower dimension"=8, "Higher Dimension"=8)
overlap_colname <- rep(c("High Overlap"=4, "Low Overlap"=4), 2)
snr_colname <- rep(c("High SNR"=2, "Low SNR"=2), 4)
yalpha_colname <- rep(1, 16)
names(yalpha_colname) <- rep(c("R", "L"), 8)

kable(rmse_table, format="latex", align="c", escape=FALSE) %>%
    add_header_above(c(" ", yalpha_colname)) %>% 
    add_header_above(c(" ", snr_colname)) %>% 
    add_header_above(c(" ", overlap_colname)) %>% 
    add_header_above(c(" ", dimension_colname))
