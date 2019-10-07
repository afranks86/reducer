library(tidyverse)
library(ggridges)
library(superheat)
library(patchwork)
library(RColorBrewer)

results_dir <- "results_test"

results_files <- dir(results_dir)

cols <- brewer.pal(5, "Set1")
cols_vec <- cols[c(1, 1, 1, 2, 2, 2, 3, 4, 5)]
lty_vec <- c("solid", "dashed", "dotted")[c(2, 3, 1, 2, 3, 1, 1, 1, 1)]

for(i in 1:length(results_files)) {

    file_name <- results_files[i]
    file_path <- paste0(results_dir, "/", file_name)

    print(file_name)

    load(file_path)
    params <- stringr::str_match(file_name,
                                 "results_n(\\d+)_p(\\d+)_coef([0-9])+_escale(-?\\d+\\.?\\d*)_mscale(-?\\d+\\.?\\d*)_yalpha(\\d+)_estpropensity(TRUE|FALSE)_([ATE]+)")
    n <- as.numeric(params[2])
    p <- as.numeric(params[3])
    coef <- as.numeric(params[4])
    escale <- as.numeric(params[5])
    mscale <- as.numeric(params[6])
    yalpha <- as.numeric(params[7])
    estimated_propensity <- as.logical(params[8]
                                       )
    estimand <- params[9]
    
    rmse_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))
    bias_mat <- apply(results_array - true_ate, c(2, 3), function(x) mean(x, na.rm=TRUE))
    var_mat <- rmse_mat^2 - bias_mat^2

    
    tib <- as_tibble(t(rmse_mat))
    tib$W <- as.numeric(colnames(rmse_mat))
    rmse_plot <- tib %>% gather(key=Type, value=RMSE, -W) %>%
        ggplot() + geom_line(aes(x=W, y=RMSE, col=Type, linetype=Type)) +
        theme_bw() + theme(legend.position="none") +
        scale_color_manual(values=cols_vec) +
        scale_linetype_manual(values=lty_vec)

    tib <- as_tibble(t(bias_mat))
    tib$W <- as.numeric(colnames(bias_mat))
    bias_plot <- tib %>% gather(key=Type, value=Bias, -W) %>%
        ggplot() + geom_line(aes(x=W, y=Bias, col=Type, linetype=Type)) +
        theme_bw() + theme(legend.position="none") + geom_hline(yintercept=0, linetype=2)  +
        scale_color_manual(values=cols_vec) +
        scale_linetype_manual(values=lty_vec)

    tib <- as_tibble(t(sqrt(var_mat)))
    tib$W <- as.numeric(colnames(bias_mat))
    sd_plot <- tib %>% gather(key=Type, value=SD, -W) %>%
        ggplot() + geom_line(aes(x=W, y=SD, col=Type, linetype=Type)) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=cols_vec) +
        scale_linetype_manual(values=lty_vec)

    if(coef == 1) {
        alpha <- c(1, -1, 0, rep(0, p-3))/sqrt(2)
        beta <- c(1, 0, 1, rep(0, p-3))/sqrt(2)
        ab_dot_prod <- t(alpha) %*% beta
    } else {
        alpha <- c(1, -1, 0, rep(0, p-3))/sqrt(2)
        beta <- c(-1, 0.75, 1, rep(0, p-3))/sqrt(3)
        ab_dot_prod <- t(alpha) %*% beta
    }
        

    if(yalpha == 0) 
        plot_title <- "Ridge Regression"
    if(yalpha == 1)
        plot_title <- "Lasso Regression"
    

    if(estimated_propensity)
        plot_title <- paste(plot_title, "Estimated Propensity Score", sep=", ")
    if(!estimated_propensity)
        plot_title <- paste(plot_title, "Known Propensity Score", sep=", ")

    ## Plot
    
    plot_name <- sprintf("results_n%i_p%i_coef%i_escale%.2f_mscale%.2f_yalpha%i_estpropensity=%s_%s",
                         n, p, coef, escale, mscale, yalpha, estimated_propensity, estimand)

    rmse_plot + bias_plot + sd_plot + plot_annotation(title = plot_title,
                                                       subtitle = sprintf("n=%i, p=%i, s_T=%.1f, s_Y=%.1f, a'b = %.2f",
                                                                          n, p, escale, mscale, ab_dot_prod))
    ggsave(filename=sprintf("figs_1/%s.pdf", plot_name), width=7, height=3)
    
}


## Make MSE Tables
library(xtable)
table_indices <- sample(length(results_files), 20)
rmse_table <- matrix(nrow=length(table_indices), ncol=8)
table_rownames <- c()
count <- 1

for(i in table_indices) {

    file_name <- results_files[i]
    file_path <- paste0("results/", file_name)

    print(file_name)

    load(file_path)
    params <- stringr::str_match(file_name,
                                 "results_n(\\d+)_p(\\d+)_coef([0-9])+_escale(-?\\d+\\.?\\d*)_mscale(-?\\d+\\.?\\d*)_yalpha(\\d+)_estpropensity(TRUE|FALSE)_([ATEC]+)")

    n <- as.numeric(params[2])
    p <- as.numeric(params[3])
    coef <- as.numeric(params[4])
    escale <- as.numeric(params[5])
    mscale <- as.numeric(params[6])
    yalpha <- as.numeric(params[7])
    estimated_propensity <- as.logical(params[8])
    estimand <- as.logical(params[9])

    rmse_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))

    rmse_table[count, ] <- c(rmse_mat[c("Naive", "Regression", "IPW", "AIPW"), 1],
                           rmse_mat[c("IPW_d", "AIPW_d"), "0"],
                           rmse_mat[c("IPW_d", "AIPW_d"), "-1"])

    table_rownames <-  c(table_rownames, sprintf("n=%i, p=%i, $s_T$=%.2f", n, p ,escale))
    count <- count + 1
}    
colnames(rmse_table) <- c("Naive", "Regression", "IPW", "AIPW", "IPW-d(0)", "AIPW-d(0)", "IPW-d(-1)", "AIPW-d(-1)")

rmse_table <- apply(rmse_table, 1, function(x) {
    x <- round(x, 2)
    min_index <- which.min(x)
    x[min_index] <- paste0("BOLD", x[min_index])
    x
}) %>% t

bold.somerows <- function(x) gsub('BOLD(.*)',paste('\\\\textbf{\\1','}'),x)

row.names(rmse_table) <- paste0(table_rownames, 1:length(table_indices))
print.xtable(xtable(rmse_table),
             sanitize.text.function=bold.somerows)




file_matches <- stringr::str_extract(results_files, "results_n1000_p1000_coef1_escale(-?\\d+\\.?\\d*)_mscale(\\d+\\.?\\d*)_yalpha[0-1]+_estpropensity(TRUE|FALSE)") 
file_matches <- results_files[!is.na(file_matches)]
file_matches <- file_matches[-c(6, 8, 10, 12, 18, 20, 22, 24)]

rmse_table <- matrix(nrow=8, ncol=length(file_matches))
table_colnames <- c()
count <- 1

escale_vec <- c()
mscale_vec <- c()
yalpha_vec <- c()


for(fn in file_matches){

    file_path <- paste0("results_1/", fn)

    load(file_path)
    params <- stringr::str_match(fn,
                                 "results_n(\\d+)_p(\\d+)_coef([0-9])+_escale(-?\\d+\\.?\\d*)_mscale(-?\\d+\\.?\\d*)_yalpha(\\d+)_estpropensity(TRUE|FALSE)")

    n <- as.numeric(params[2])
    p <- as.numeric(params[3])
    coef <- as.numeric(params[4])
    escale <- as.numeric(params[5])
    mscale <- as.numeric(params[6])
    yalpha <- as.numeric(params[7])
    estimated_propensity <- as.logical(params[8])

    rmse_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))

    rmse_table[, count] <- c(rmse_mat[c("Naive", "Regression", "IPW", "AIPW"), 1],
                             rmse_mat[c("IPW_d", "AIPW_d"), "0"],
                             rmse_mat[c("IPW_d", "AIPW_d"), "-1"])
    
    table_colnames <-  c(table_colnames, paste(coef, escale, mscale, yalpha, estimated_propensity, sep=","))

    escale_vec <- c(escale_vec, escale)
    mscale_vec <- c(mscale_vec, mscale)
    yalpha_vec <- c(yalpha_vec, yalpha)
    
    count <- count + 1
}

rownames(rmse_table) <- c("Naive", "Regression", "IPW", "AIPW", "IPW-d(0)", "AIPW-d(0)", "IPW-d(-1)", "AIPW-d(-1)")

rmse_table <- apply(rmse_table, 2, function(x) {
    x <- round(x, 2)
    min_index <- which.min(x)
    x[min_index] <- paste0("BOLD", x[min_index])
    x
})

addtorow <- list()
addtorow$pos <- list(0, 1)
addtorow$command <- c(paste0(paste0('& \\multicolumn{12}{cc}{', unique(escale_vec), '}', collapse=''), '\\\\'),
                      paste0(paste0('& \\multicolumn{12}{c}{', 1:2, '}', collapse=''), '\\\\'))

colnames(rmse_table) <- table_colnames
print.xtable(xtable(rmse_table),
             sanitize.text.function=bold.somerows)

colnames(rmse_table) <- table_colnames
print.xtable(xtable(rmse_table),
             add.to.row = addtorow,
             sanitize.text.function=bold.somerows,
             include.colnames = FALSE)

paste0(escale_vec, collapse= " & ")
paste0(mscale_vec, collapse= " & ")
paste0(ifelse(yalpha_vec==0, "R", "L"), collapse=" & ")
paste0(rep(c("F", "T"), 8), collapse=" & ")
