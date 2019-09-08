library(tidyverse)
library(ggridges)
library(superheat)
library(patchwork)
as_tibble(results_array[, , "0.6"]) %>%
    gather(key=Type, value=Estimate) %>%
    ggplot() +
    geom_density_ridges(aes(x=Estimate, y=Type, fill=Type), stat="binline") +
    geom_vline(xintercept=5) + theme_bw()


rmse_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))
bias_mat <- apply(results_array - true_ate, c(2, 3), function(x) mean(x, na.rm=TRUE))
var_mat <- rmse_mat^2 - bias_mat^2


tib <- as_tibble(t(rmse_mat))
tib$value <- as.numeric(colnames(rmse_mat))
rmse_plot <- tib %>% gather(key=Type, value=RMSE, -value) %>%
    ggplot() + geom_line(aes(x=value, y=RMSE, col=Type)) +
    theme_bw() + theme(legend.position="none")

tib <- as_tibble(t(bias_mat))
tib$value <- as.numeric(colnames(bias_mat))
bias_plot <- tib %>% gather(key=Type, value=Bias, -value) %>%
    ggplot() + geom_line(aes(x=value, y=Bias, col=Type)) +
    theme_bw() + theme(legend.position="none")

tib <- as_tibble(t(var_mat))
tib$value <- as.numeric(colnames(bias_mat))
var_plot <- tib %>% gather(key=Type, value=Var, -value) %>%
    ggplot() + geom_line(aes(x=value, y=Var, col=Type)) + theme_bw() 


rmse_plot + bias_plot + var_plot

cor(results_array[, "IPW_d",  ])
superheat::superheat(cor(results_array[, "AIPW_d",  ]))

cor(results_array[, ,  "0.9"])
cor(results_array[, ,  "0.7"])