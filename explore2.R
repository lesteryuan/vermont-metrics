## 12.4.2023
## initial exploration

## edit bcnt file and generate otu file
explore <- function(bug.data) {

    ## loaded CSV files into R as bug.data and site.data.

    ## make sure sample id is a factor
    bug.data$bio_sample <- factor(bug.data$bio_sample)

    ## replace text "N/A" with R NA
    fnames <- c("Order", "Family", "SubFamilyOrTribe", "GenusGroup",
                "Genus","SpeciesGroup", "Species")
    for (i in fnames) {
        incvec <- bug.data[, i] == "N/A"
        bug.data[incvec,i] <- NA
    }

    ## change "sp" in species field to NA
    incvec <- tolower(bug.data$Species) == "sp" |
        bug.data$Species == "sp a" | bug.data$Species == "sp b" |
            bug.data$Species == "" | bug.data$Species == "w/o setae" |
                  bug.data$Species == "spa" |
                    bug.data$Species == "uid" |
                        bug.data$Species == "imm" |
                            bug.data$Species == "uid wsetae" |
                                tolower(bug.data$Species) == "group"
    bug.data$Species[incvec] <- NA


#    incvec <- bug.data$SpeciesGroup == "uid" |
#        tolower(bug.data$SpeciesGroup) == "group"
#    print(table(bug.data$SpeciesGroup))
#    bug.data$SpeciesGroup[incvec] <- NA

    ## set ambiguous Genus to NA
    incvec <- regexpr("genus", tolower(bug.data$Genus)) != -1
    bug.data$Genus[incvec] <- NA
    incvec <- bug.data$Genus == "UID" | bug.data$Genus == "UNID" |
        tolower(bug.data$Genus) == "group"
    incvec[is.na(incvec)] <- F
    bug.data$Genus[incvec] <- NA

    incvec <- tolower(bug.data$GenusGroup) == "group"
    incvec[is.na(incvec)] <- FALSE
    bug.data$GenusGroup[incvec] <- NA

    ## Combine genus, genusgroup, subfamily/tribe into one taxon
    ## such that taxa not identified to genus are recorded as genusgroup
    ## and so on.
    taxon <- bug.data$Genus

    incvec <- !is.na(bug.data$Species) & !is.na(bug.data$Genus)
    taxon[incvec] <- paste(bug.data$Genus[incvec], bug.data$Species[incvec])

    incvec <- is.na(bug.data$Species) & ! is.na(bug.data$SpeciesGroup) &
        ! is.na(bug.data$Genus)
    taxon[incvec] <- paste(bug.data$Genus[incvec], bug.data$SpeciesGroup[incvec])

    print(sum(is.na(taxon)))
    incvec <- is.na(taxon)
    taxon[incvec] <- bug.data$GenusGroup[incvec]

    print(sum(is.na(taxon)))
    incvec <- is.na(taxon)
    taxon[incvec] <- bug.data$SubFamilyOrTribe[incvec]

    print(sum(is.na(taxon)))
    incvec <- is.na(taxon)
    taxon[incvec] <- bug.data$Family[incvec]

    print(sum(is.na(taxon)))
    bug.data$taxon <- taxon

    incvec <- is.na(bug.data$taxon)
    bug.data <- bug.data[!incvec,]

    incvec <- bug.data$density_m2 > 0
    ## make zero densities a small number to
    ## account for a presence
    bug.data$density_m2[incvec] <- 0.01

    require(bio.infer)
    bcnt.tax <- get.taxonomic(bug.data[, c("bio_sample", "taxon", "density_m2")])

    source("get.otuloc.R")

    bcnt.otu <- get.otuloc(bcnt.tax, outputFile = T)

    ## further edits to sum.otu.txt are required externally
    ## to pick out a few more species rather than genus ids
    ## these are loaded back up with load.revised.otu
    return(bcnt.otu)
}


## make combined data frame with site-species and envdata
makefinaldf <- function(ss, site.data) {

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

    ## select taxa that occurred in at 20 samples for further analysis
    numocc <- apply(ss[, tnames], 2, function(x) sum(x>0))
    numocc.sav <- numocc[numocc > 20]
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

    ## turbidity has a bunch of measurements that are zero
    ## set these to half the minimum positive value that was detected
    ## to allow for log transform
    incvec <- ss$turbidity == 0
    incvec[is.na(incvec)] <- F
    minval <- 0.5*min(ss$turbidity[!incvec], na.rm = T)
    ss$turbidity[incvec] <- minval
    ss$turbidity <- log(ss$turbidity)
    attr(ss, "taxa") <- tnames  ## assign selected taxa as attributes
    return(ss)
}

runwa <- function(ss, varname) {
    
    tnames <- attr(ss, "taxa")

    pred.wa <- matrix(NA, ncol = length(varname), nrow = nrow(ss))
    dimnames(pred.wa)[[2]] <- varname
    #        png(width = 6, height = 4, units = "in", res = 600, pointsize = 8,
#            file = "taxop.png")
    dev.new()
    par(mar = c(4,4,3,1), mfrow = c(2,3), mgp = c(2.3,1,0))
    for (j in 1:length(varname)) {
        incvec <- ! is.na(ss[, varname[j]])
        print(sum(incvec))
        ss0 <- ss[incvec,]
        
        opt <- rep(NA, times = length(tnames))
        names(opt) <- tnames
        for (i in 1:length(tnames)) {
            selvec <- ss0[, tnames[i]] > 0
            opt[i] <- mean(ss0[selvec, varname[j]], na.rm = T)
            
        }
        inf <- rep(NA, times = nrow(ss0))
        for (i in 1:nrow(ss0)) {
            selvec <- ss0[i, tnames] > 0
            inf[i] <- mean(opt[selvec])
        }
        pred.wa[incvec,j] <- inf
        modfit <- lm(ss0[, varname[j]]~inf)
        
        plot(inf, ss0[, varname[j]], pch = 21,
             col = "grey", bg = "white",
             xlab = paste(lab0[j], "(Predicted)"),
             ylab = paste(lab0[j], "(Observed)"), axes = F)
        mtext(paste("R2 =", round(summary(modfit)$r.squared, digits = 2)),
              side = 3, line = 0,
              cex = 0.8)
        if (logt[j]) {
            logtick.exp(0.001, 10, c(1,2), c(T,F))
        }
        else {
            axis(1)
            axis(2)
            box(bty = "l")
        }
        abline(modfit)
    }
    return(pred.wa)
}

## run RF model for inference
runRF <- function(ss, varname) {
    require(ranger)
    tnames <- attr(ss, "taxa")

    ## set up storage locations
    predsav <- matrix(NA, ncol = length(varname), nrow = nrow(ss))
    dimnames(predsav)[[2]] <- varname
    pefflist <- as.list(rep(NA, times = length(varname)))
    names(pefflist) <- varname
                        
    ## fit models
    for (j in 1:length(varname)) {

        ## drop missing measurements
        incvec <- ! is.na(ss[, varname[j]])
        ss0 <- ss[incvec,]

        domtry <- F ## optimize selection of mtry
        if (domtry) {
            kmtry <- seq(1, 1.6, by = 0.1)
            ## mtry = 21 taxa works well for turbidity
            err <- rep(NA, times = length(kmtry))
            for (i in 1:length(kmtry)) {
                mod <- ranger(data = ss0[, c(varname[j], tnames)],
                              dependent.variable.name = varname[j],
                              num.trees = 6000,mtry = round(sqrt(length(tnames))*kmtry[i]),
                              importance = "permutation")
                err[i] <- mod$prediction.error
            }
            plot(kmtry, err)
            stop()
        }
        
        set.seed(10)
        ## fit RF 20 times to get multiple assessments of significant taxa
        nrun <- 20
        namesav <- as.list(rep(NA, times =nrun))
        for (i in 1:nrun) {
            mod.imp <- ranger(data = ss0[, c(varname[j], tnames)],
                              dependent.variable.name = varname[j],
                              num.trees = 5000,
                              importance = "impurity_corrected")
#        print(mod.imp)
        ## select statistically significant taxa (initial list of candidates)
            imp0 <- importance_pvalues(mod.imp, method = "janitza")
            namesav[[i]] <- dimnames(imp0)[[1]][imp0[,2] < 0.05]
        }
        namesall <- unlist(namesav)

        ## pick taxa that are significant in at least half the runs
        tout <- table(namesall)
        namesav <- names(tout)[tout > 0.5*nrun]
        print(length(namesav))

        ## rerun model with selected taxa
        mod <- ranger(data = ss0[, c(varname[j], namesav)],
                      dependent.variable.name = varname[j], num.trees = 5000,
                      importance = "permutation")
        print(mod)

        predsav[incvec,j] <- mod$predictions

        ## calculate opt for each taxon
        require(pdp)
        peff <- rep(NA, times = length(namesav))
        names(peff) <- namesav
        for (i in 1:length(namesav)) {
            cat(i, " ")
            if (floor(i/10) == i/10) cat("\n")
            flush.console()
            pout <- partial(mod, pred.var = namesav[i], plot = FALSE)
            peff[i] <- (pout$yhat[2] - pout$yhat[1])/
                diff(range(predsav[,j], na.rm = T))
        }
        cat("\n")

        pefflist[[j]] <- peff

    }

    return(list(predsav, pefflist))
}

postp <- function(ss, varname,  pred.RF, lab0, logt) {
    dev.new()
    par(mar = c(4,4,3,1), mfrow = c(2,3), mgp= c(2.3,1,0))
    for (j in 1:length(varname)) {
        plot(pred.RF[[1]][, varname[j]], ss[, varname[j]], pch = 21,
             col = "grey", bg = "white",
             xlab = paste(lab0[j], "(Predicted)"),
             ylab = paste(lab0[j], "(Observed)"), axes = F)
        mod <- lm(ss[, varname[j]] ~ pred.RF[[1]][, varname[j]])
        mtext(paste("R2 =", round(summary(mod)$r.squared, digits = 2)),
              side = 3, line = 0,
              cex = 0.8)

        if (logt[j]) {
            logtick.exp(0.001, 10, c(1,2), c(F,F))
        }
        else {
            axis(1)
            axis(2)
            box(bty = "l")
        }
        abline(0,1, lty = "dashed")
        abline(mod)
    }

    dev.new()
    par(mar = c(4,4,3,1), mfrow = c(2,3), mgp= c(2.3,1,0))    
    ## simple inference using RF optima
    for (j in 1:length(varname)) {
        peff <- pred.RF[[2]][[j]]
        infout <- as.matrix(ss[, names(peff)]) %*% peff
        plot(infout, ss[, varname[j]], pch = 21,
             col = "grey", bg = "white",
             xlab = paste(lab0[j], "(Predicted)"),
             ylab = paste(lab0[j], "(Observed)"), axes = F)
        mod <- lm(ss[, varname[j]] ~ infout)
        mtext(paste("R2 =", round(summary(mod)$r.squared, digits = 2)),
              side = 3, line = 0,
              cex = 0.8)

        if (logt[j]) {
            logtick.exp(0.001, 10, c(1,2), c(F,F))
        }
        else {
            axis(1)
            axis(2)
            box(bty = "l")
        }
        abline(0,1, lty = "dashed")
        abline(mod)
    }
    stop()

    plot(peff, opt[names(peff)])
    stop()
}


## cluster predictor variables
cluster.env <- function(df1) {

    names0 <- c("latitude", "longitude", "area_km2", "elev_ft",
                "conductivity", "chloride_mg", "tp_ug", "nitrate_mg",
                "alkalinity", "sulfate_mg", "tal_ug", "tfe_ug",
                "turbidity",
                "watershed_imperv", "watershed_developed",
                "riparian_developed", "watershed_agric",
                "riparian_agric", "watershed_forest",
                "riparian_forest", "watershed_wetland",
                "riparian_wetland", "riparian_natural",
                "finesed_percent", "gravel_percent",
                "embedded_percent")
    df1 <- na.omit(df1[, names0])

    cor0 <- 1-abs(cor(df1[, names0], method = "spearman"))

    require(cluster)
    clust0 <- agnes(cor0, diss = T, method = "average")

    lab0 <- names0

    png(width = 4, height = 4, pointsize = 7, units = "in", res = 600,
            file = "dend.png")
    par(mar = c(2,1,1,4), mgp = c(2.3,1,0))
    plot(clust0,labels = lab0,  which.plots = 2,
         main = "", ylab = "" , xlab = "", sub = "", axes
         = F)
    axis(4, at = seq(0, 1, by = 0.2), lab = seq(1, 0, by = -0.2))
    mtext(expression(abs(italic(r[s]))),side = 4, line = 2.3)

    dev.off()
    return()
}


#bcnt.otu <- explore(bug.data1229)

## make sure to run these at species (from Aaron)
## optioservus, stenelmis, polypedilum, rheotanytarsus, potthastia,
## eukiefferiella, tvetenia, orthocladius/euorthocladius, simulium,
## baetis, ephemerella*, maccaffertium*, brachycentrus, hydropsyche,
## chimarra, rhyacophila, acroneuria, paragnetina, pteronarcys,

## *keeping all macaffertium rather than just m. vicarium yields a
## better model for salinity
## ephemerella is an indicator for salinity, but not when only species
## level data are used.

#cluster.env(site.data)

#ssall <- makefinaldf(ss, site.data)

varname <- c("turbidity", "silt_rating", "alkalinity", "chloride",
             "tp_ug")
lab0 <- c("Turbidity", "Silt rating", "Alkalinity", "Chloride",
          "Total P")
logt <- c(T, F, T, T, T)

#pred.wa <-runwa(ssall, varname)

#pred.RF <- runRF(ssall, varname)

postp(ssall, varname, pred.RF, lab0, logt)
