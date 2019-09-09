library(ftidyverse)
library(ggridges)
library(superheat)
library(patchwork)
results_files <- dir("results")

for(i in 1:length(results_files)) {

    file_name <- results_files[i]
    file_path <- paste0("results/", file_name)

    print(file_name)

    load(file_path)
    params <- stringr::str_match(file_name,
                                 "results_n(\\d+)_p(\\d+)_coef([0-9])+_escale(-?\\d+\\.?\\d*)_mscale(-?\\d+\\.?\\d*)_yalpha(\\d+)_estpropensity(TRUE|FALSE)")

    n <- as.numeric(params[2])
    p <- as.numeric(params[3])
    coef <- as.numeric(params[4])
    escale <- as.numeric(params[5])
    mscale <- as.numeric(params[6])
    yalpha <- as.numeric(params[7])
    estimated_propensity <- as.logical(params[8]
)
    rmse_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))
    bias_mat <- apply(results_array - true_ate, c(2, 3), function(x) mean(x, na.rm=TRUE))
    var_mat <- rmse_mat^2 - bias_mat^2


    tib <- as_tibble(t(rmse_mat))
    tib$W <- as.numeric(colnames(rmse_mat))
    rmse_plot <- tib %>% gather(key=Type, value=RMSE, -W) %>%
        ggplot() + geom_line(aes(x=W, y=RMSE, col=Type)) +
        theme_bw() + theme(legend.position="none")

    tib <- as_tibble(t(bias_mat))
    tib$W <- as.numeric(colnames(bias_mat))
    bias_plot <- tib %>% gather(key=Type, value=Bias, -W) %>%
        ggplot() + geom_line(aes(x=W, y=Bias, col=Type)) +
        theme_bw() + theme(legend.position="none") + geom_hline(yintercept=0, linetype=2)

    tib <- as_tibble(t(sqrt(var_mat)))
    tib$W <- as.numeric(colnames(bias_mat))
    sd_plot <- tib %>% gather(key=Type, value=SD, -W) %>%
        ggplot() + geom_line(aes(x=W, y=SD, col=Type)) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))

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
    
    plot_name <- sprintf("results_n%i_p%i_coef%i_escale%.2f_mscale%.2f_yalpha%i_estpropensity=%s",
                         n, p, coef, escale, mscale, yalpha, estimated_propensity)

    rmse_plot + bias_plot + sd_plot + plot_annotation(title = plot_title,
                                                       subtitle = sprintf("n=%i, p=%i, escale=%.1f, mscale=%.1f, a'b = %.2f",
                                                                          n, p, escale, mscale, ab_dot_prod))
    ggsave(filename=sprintf("figs/%s.pdf", plot_name), width=14)



    
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
                                 "results_n(\\d+)_p(\\d+)_coef([0-9])+_escale(-?\\d+\\.?\\d*)_mscale(-?\\d+\\.?\\d*)_yalpha(\\d+)_estpropensity(TRUE|FALSE)")

    n <- as.numeric(params[2])
    p <- as.numeric(params[3])
    coef <- as.numeric(params[4])
    escale <- as.numeric(params[5])
    mscale <- as.numeric(params[6])
    yalpha <- as.numeric(params[7])
    estimated_propensity <- as.logical(params[8])

    rmse_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))

    rmse_table[count, ] <- c(rmse_mat[c("Naive", "Regression", "IPW", "AIPW"), 1],
                           rmse_mat[c("IPW_d", "AIPW_d"), "0"],
                           rmse_mat[c("IPW_d", "AIPW_d"), "-1"])

    table_rownames <-  c(table_rownames, sprintf("n=%i, p=%i, escale=%.2f", n, p ,escale))
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

