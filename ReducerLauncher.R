library(lubridate)
RSCRIPT_ALIAS <- "Rscript"

iters <- 10

np <- list(c(100, 100), c(200, 100), c(1000, 500), c(1000, 1000))
coef_settings <- c(1, 2)
mscale <- c(-10, 5, 5, 10)
escale <- c(1, 3)
y_alpha <- c(0, 1)
EST_PROPENSITY <- c(FALSE, TRUE)

option_names <- c('n', 'p',   'coef', 'escale', 'mscale', 'y_alpha', 'est_propensity')
option_types <- c('%i', '%i', '%i',   '%.1f',   '%.1f',   '%i',      '%s')
option_fstring <- paste('--', option_names, '=', option_types, collapse=' ', sep='')

script_fstring <- paste(RSCRIPT_ALIAS, "ReducerSims.R", option_fstring, sprintf("--iters=%i", iters))

logfile_fstring <- "logs/log_n%i_p%i_coef%i_escale%.1f_mscale%.1f_yalpha%i_estpropensity%s_%s.log"

for(aspects in np){
    for(cc in coef_settings){
        for(ee in escale){
            for(mm in mscale){
                for(aa in y_alpha){
                    for(est in EST_PROPENSITY){
                        call <- sprintf(script_fstring, aspects[1], aspects[2], cc, ee, mm, aa, est)
                        print(call)
                        logfile <- sprintf(logfile_fstring, aspects[1], aspects[2], cc, ee, mm, aa, est,
                                           now())
                        system(paste(call, ">", logfile, "2>&1"), wait=FALSE)
                    }
                }
            }
        }
    }
}