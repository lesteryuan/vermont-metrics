## 1.31.2024
## try multivariable RF to alleviate correlation

explore2 <- function(ss, site.data, keytaxa = NULL) {

    ## replace slashes with underscores for later formulas
    tnames <- names(ss)[-1]
    tnames.rev <- gsub("/", "_", tnames)
    print(tnames.rev)
    names(ss) <- c(names(ss)[1], tnames.rev)
    tnames <- tnames.rev

    ## make ss into presence/absence
    for (i in tnames) {
        incvec <- ss[, i] > 0
        ss[incvec,i] <- 1
    }

    site.data$bio_sample <- factor(site.data$bio_sample)
    ## drop 1 duplicate bio_sample
    site.data <- site.data[unique(site.data$bio_sample),]
    ## select kicknet
    incvec <- site.data$sample_type == "KN"
    site.data <- site.data[incvec,]

    ss <- merge(ss, site.data,by = "bio_sample")
    print(dim(ss))

    ## select taxa that occurred in at 30 samples for further analysis
    numocc <- apply(ss[, tnames], 2, function(x) sum(x>0))
    numocc.sav <- numocc[numocc > 30]
    tnames <- names(numocc.sav)
    print(tnames)

    ## some log-transforms of the chemistry
    ss$conductivity <- log(ss$conductivity)
    ss$tp_ug <- log(ss$tp_ug)
    ss$nitrate_mg <- log(ss$nitrate_mg) # detection limit issues
    ss$chloride <- log(ss$chloride_mg) # detection limit issues
    ## drop one outlier alkalinity
    incvec <- ss$alkalinity > 0.1
    incvec[is.na(incvec)] <- T
    ss <- ss[incvec,]
    ss$alkalinity <- log(ss$alkalinity)

 #   plot(ss$alkalinity, ss$chloride)

    ## turbidity has a bunch of measurements that are zero
    ## set these to half the minimum positive value that was detected
    ## to allow for log transform
    incvec <- ss$turbidity == 0
    incvec[is.na(incvec)] <- F
    minval <- 0.5*min(ss$turbidity[!incvec], na.rm = T)
    ss$turbidity[incvec] <- minval
    ss$turbidity <- log(ss$turbidity)

    require(randomForestSRC)
    varname <- c("turbidity", "silt_rating", "alkalinity", "chloride",
                 "tp_ug")
    lab0 <- c("Turbidity", "Silt rating", "Alkalinity", "Chloride",
              "Total P")
    logt <- c(T, F, T, T, T)

    ## pairs plot to examine covariance
#    dev.new()
#    pairs(ss[, varname])
#    print(cov(ss[, varname], use ="pair"))

    ## set up storage locations
    predsav <- matrix(NA, ncol = length(varname), nrow = nrow(ss))
    dimnames(predsav)[[2]] <- varname
    peff <- matrix(NA, ncol = length(varname), nrow = length(tnames))
    dimnames(peff)[[1]] <- tnames
    dimnames(peff)[[2]] <- varname
    imp <- matrix(NA, ncol = length(varname), nrow = length(tnames))
    dimnames(imp)[[1]] <- tnames
    dimnames(imp)[[2]] <- varname

    
    incvec <- ! is.na(ss[, "chloride"]) & ! is.na(ss[, "tp_ug"])
    print(sum(incvec))
    ss0 <- ss[incvec,]

    ## scale variables
    ss0$chloride.sc <- (ss0$chloride - mean(ss0$chloride))/sd(ss0$chloride)
    ss0$tp_ug.sc <- (ss0$tp_ug - mean(ss0$tp_ug))/sd(ss0$tp_ug)
    
    mod <- rfsrc(Multivar(chloride.sc, tp_ug.sc) ~ .,
                 data = ss0[, c("chloride.sc", "tp_ug.sc", tnames)], 
                 splitrule = "mahalanobis")

    pout <- get.mv.predicted(mod)
    print(str(pout))

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2))
    plot(pout[,1], ss0[, "chloride.sc"])
    abline(0,1)
    print(summary(lm(ss0[, "chloride.sc"] ~ pout[,1])))
    plot(pout[,2], ss0[, "tp_ug.sc"])
    abline(0,1)
    print(summary(lm(ss0[, "tp_ug.sc"] ~ pout[,2])))

    print(cor(pout[,1], pout[,2]))

    ## multivariable inferences are more correlated than single variable
    ## models...
    stop()
}


explore2(ss, site.data, keytaxa)
