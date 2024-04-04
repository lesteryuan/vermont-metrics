## look for interactions between taxa

f_interact <- function(ss, var, tlist, dfout = NULL) {
    require(rpart)

    ## run rpart for every pairwise combo
    names0 <- names(tlist)

    if (is.null(dfout)) {
        ic <- 0
        ncomb <- (length(names0)^2 - length(names0))/2
        delt <- rep(NA, times = ncomb)
        isav <- rep(NA, times = ncomb)
        jsav <- rep(NA, times = ncomb)
        for (i in 1:(length(names0)-1)) {
            for (j in (i+1):length(names0)) {
                ic <- ic + 1
                formstr <- paste(var, "~", names0[i], "+", names0[j])
                mod <- rpart(as.formula(formstr), ss, control = list(minbucket = 5))

                new.data <- data.frame(c(1,1,0,0), c(1,0, 1,0))
                names(new.data) <- c(names0[i], names0[j])
                predout <- predict(mod, new.data)
                delt[ic] <- predout[4] - predout[3] - (predout[2] - predout[1])
                isav[ic] <- i
                jsav[ic] <- j
            }
        }
        dfout <- data.frame(i0 = isav, j0 = jsav, delt = delt)
        return(dfout)
    }

    else {
        delt.sc <- abs(dfout$delt/diff(quantile(ss[, var], prob = c(0.05, 0.95),
                                            na.rm = T)))
        incvec <- delt.sc > 0.5
        print(sum(incvec))
        tfreq <- rev(sort(table(c(dfout$i0[incvec], dfout$j0[incvec]))))
        print(names0[as.numeric(names(tfreq)[1:3])])
    }
    
}
#dfout <- f_interact(ssall, "watershed_imperv", pred.RF.all[[2]][["watershed_imperv"]])
#dfout.turb <- f_interact(ssall, "turbidity", pred.RF.all[[2]][["turbidity"]])
f_interact(ssall, "watershed_imperv", pred.RF.all[[2]][["watershed_imperv"]],
           dfout = dfout)
