library(lubridate)
library(parallel)
RSCRIPT_ALIAS <- "/opt/R/3.5.3/bin/Rscript"

iters <- 500 

#np <- list(c(100, 100), c(200, 100), c(1000, 500), c(1000, 1000))
#coef_settings <- c(1, 2)
#escale <- c(1, 4)
#mscale <- c(-10, 5, 5, 10)
#y_alpha <- c(0, 1)
#EST_PROPENSITY <- c(FALSE, TRUE)

np <- list(c(500, 100), c(500, 1000))
coef_settings <- c(1, 3) ## sparsity of outcome model
ab_dp <- c(0.75) ## Dot product between alpha and beta
escale <- c(1, 4)
mscale <- c(2, 5)
y_alpha <- c(0, 1)
EST_PROPENSITY <- c(FALSE)

all_settings <- expand.grid(np, coef_settings, escale, mscale, y_alpha, EST_PROPENSITY)

option_names <- c('n', 'p',   'coef', 'escale', 'mscale', 'y_alpha', 'est_propensity')
option_types <- c('%i', '%i', '%i',   '%.1f',   '%.1f',   '%i',      '%s')
option_fstring <- paste('--', option_names, '=', option_types, collapse=' ', sep='')

script_fstring <- paste(RSCRIPT_ALIAS, "ReducerSims.R", option_fstring, sprintf("--iters=%i", iters))

logfile_options <- paste(option_names, option_types, collapse='_', sep='')
logfile_fstring <- paste("logs/log_", logfile_options, "_%s.log", sep='')

run_setting <- function(row){
    row <- unlist(row)
    n <- row[1]
    p <- row[2]
    cc <- row[3]
    ee <- row[4]
    mm <- row[5]
    aa <- row[6]
    est <- as.logical(row[7])
    
    call <- sprintf(script_fstring, n, p, cc, ee, mm, aa, est)
    print(call)
    logfile <- sprintf(logfile_fstring, n, p, cc, ee, mm, aa, est,
                       gsub(" ", "", now(), fixed=TRUE))
    system(paste(call, ">", logfile, "2>&1"))
}


retcodes <- mclapply(1:nrow(all_settings),
         function(i){
             run_setting(all_settings[i,])
         }, mc.cores=detectCores())

print(retcodes)
save(retcodes, all_settings,
     file=paste("logs/experiment_exit_status_",
                gsub(" ", "", now(), fixed=TRUE), ".log", sep=""))
