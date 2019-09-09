library(lubridate)
library(parallel)
RSCRIPT_ALIAS <- "/opt/R/3.5.3/bin/Rscript"

iters <- 100 

np <- list(c(100, 100), c(200, 100), c(1000, 500), c(1000, 1000))
coef_settings <- c(1, 2)
escale <- c(1, 4)
mscale <- c(-10, 5, 5, 10)
y_alpha <- c(0, 1)
EST_PROPENSITY <- c(FALSE, TRUE)

all_settings <- expand.grid(np, coef_settings, escale, mscale, y_alpha, EST_PROPENSITY)

option_names <- c('n', 'p',   'coef', 'escale', 'mscale', 'y_alpha', 'est_propensity')
option_types <- c('%i', '%i', '%i',   '%.1f',   '%.1f',   '%i',      '%s')
option_fstring <- paste('--', option_names, '=', option_types, collapse=' ', sep='')

script_fstring <- paste(RSCRIPT_ALIAS, "ReducerSims.R", option_fstring, sprintf("--iters=%i", iters))

logfile_fstring <- "logs/log_n%i_p%i_coef%i_escale%.1f_mscale%.1f_yalpha%i_estpropensity%s_%s.log"

retcodes <- mclapply(1:nrow(all_settings),
         function(i){
             n <- all_settings[i,1][[1]][1]
             p <- all_settings[i,1][[1]][2]
             cc <- all_settings[i,2]
             ee <- all_settings[i,3]
             mm <- all_settings[i,4]
             aa <- all_settings[i,5]
             est <- all_settings[i,6]
             
             call <- sprintf(script_fstring, n, p, cc, ee, mm, aa, est)
             print(call)
             logfile <- sprintf(logfile_fstring, n, p, cc, ee, mm, aa, est,
                                gsub(" ", "", now(), fixed=TRUE))
             system(paste(call, ">", logfile, "2>&1"))
         }, mc.cores=detectCores())

print(retcodes)
