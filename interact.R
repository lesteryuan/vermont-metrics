## look for interactions between taxa

## rerun RF model without the strong interactors and see if the
## r2/r20 gets closer to 1.
## will have to speculate on why interactions are more important
## for impervious surfaces compared to chemical stressors.

f_interact <- function(ss, var, tlist, dfout = NULL, pred.RF = NULL) {
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
                mod <- rpart(as.formula(formstr), ss, control = list(minbucket = 20))

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
        delt.sc <- abs(dfout$delt/diff(range(ss[, var],na.rm = T)))

        dev.new()
        hist(delt.sc)
        print(summary(delt.sc))

        incvec <- delt.sc > 0.25
        print(sum(incvec))
        tfreq <- rev(sort(table(c(dfout$i0[incvec], dfout$j0[incvec]))))
        print(tfreq)
        ip <- as.numeric(names(tfreq)[1:3])
        print(names0[ip])

        dfout.red <- dfout[incvec,]

        require(maps)
        require(mapproj)
        dev.new()
        map("state", proj = "albers", region = "vermont", par = c(30,40))
        pout <- mapproject(ss$longitude, ss$latitude, proj = "")
        incvec <- ss$CALOPTERYX == 1
#        incvec2 <- ss$CALOPTERYX == 1 & ss$EPHEMERELLA == 0
        incvec2 <- ss$CALOPTERYX == 1 & ss$OULIMNIUS == 0
        print(sum(incvec2))
        points(pout$x[incvec], pout$y[incvec])
        points(pout$x[incvec2], pout$y[incvec2], pch = 16)

        print(ss[incvec2, c("latitude", "longitude")])

        for (i in ip[1]) {
            incvec<- dfout.red$i0 ==i | dfout.red$j0 == i
            dftemp <- dfout.red[incvec,]

            par(mar = c(4,4,1,1), mfrow = c(2,3))
            icount <- 0
            for (j in 1:nrow(dftemp)) {
                icount <- icount + 1
                formstr <- paste(var, "~", names0[dftemp$i0[j]], "+",
                                 names0[dftemp$j0[j]])
                mod <- rpart(as.formula(formstr), ss,
                             control = list(minbucket = 5))
                plot(mod, uniform = T, compress = T, margin = 0.3)
                text(mod)


                ifac <- interaction(ss[, names0[dftemp$i0[j]]],
                                    ss[, names0[dftemp$j0[j]]])
                boxplot(split(ss[, var], ifac))
                stop()

                if (round(icount/6) == icount/6) {
                    dev.new()
                    par(mar = c(4,4,1,1), mfrow = c(2,3))
                }
            }
            stop()
        }


    }

}
#dfout <- f_interact(ssall, "watershed_imperv", pred.RF.all[[2]][["watershed_imperv"]])
#dfout.turb <- f_interact(ssall, "turbidity", pred.RF.all[[2]][["turbidity"]])
f_interact(ssall, "watershed_imperv", pred.RF.all[[2]][["watershed_imperv"]],
           dfout = dfout, pred.RF = pred.RF.all)
