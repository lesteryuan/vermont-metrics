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


## start analysis from ss from bcnt.otu
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

#    dev.new()
#    png(width = 5, height = 2.5, pointsize = 7, units = "in", res = 600,
#        file = "impplot.png")
#    par(mar = c(4,4,2,2),las=1, mfrow = c(1,2), mgp = c(2.3,1,0))
#    png(width = 6, height = 4, units = "in", res = 600, pointsize = 8,
#        file = "taxop.png")
#    par(mar = c(4,4,3,1), mfrow = c(2,3), mgp = c(2.3,1,0))
#    png(width = 5, height = 5, units = "in", res = 600, pointsize = 8,
#        file = "taxop.png")
#    par(mar = c(4,10,1,1), mfrow = c(1,2), mgp = c(2.3,1,0), las = 1)

    require(ranger)
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

    print(cor(ss[, varname], use = "pair"))
    stop()

    ## run weighted average to compare covariance
    ## weighted averages are super-correlated, so RF is a big
    ## improvement
    dowa <- F
    if (dowa) {
        pred.wa <- matrix(NA, ncol = length(varname), nrow = nrow(ss))
        dimnames(pred.wa)[[2]] <- varname
        png(width = 6, height = 4, units = "in", res = 600, pointsize = 8,
            file = "taxop.png")
        par(mar = c(4,4,3,1), mfrow = c(2,3), mgp = c(2.3,1,0))
        for (j in 1:length(varname)) {
            opt <- rep(NA, times = length(tnames))
            names(opt) <- tnames
            for (i in 1:length(tnames)) {
                selvec <- ss[, tnames[i]] > 0
                opt[i] <- mean(ss[selvec, varname[j]], na.rm = T)

            }
            inf <- rep(NA, times = nrow(ss))
            for (i in 1:nrow(ss)) {
                selvec <- ss[i, tnames] > 0
                inf[i] <- mean(opt[selvec])
            }
            pred.wa[,j] <- inf
            modfit <- lm(ss[, varname[j]]~inf)

            plot(inf, ss[, varname[j]], pch = 21,
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
        dev.off()
        save(pred.wa, file = "pred.wa.rda")
        stop()
    }

    ## these are prediction error cutoffs that define which
    ## taxa to retain in the model. Identified by trying different
    ## prediction error cutoffs to minimize overall error
    cutsav <- c(0.006, 0.015, 0.007, 0.015, 0.003)
    names(cutsav) <- varname

    ## keytaxa are those selected by cutsav
    if (is.null(keytaxa)) {
        keytaxa <- as.list(rep(NA,times = length(varname)))
    }

    ## fit models
    for (j in 1:length(varname)) {
        incvec <- ! is.na(ss[, varname[j]])
        print(sum(incvec))
        ss0 <- ss[incvec,]

        set.seed(10)
        ## fit initial RF to predict varname
        mod.imp <- ranger(data = ss0[, c(varname[j], tnames)],
                      dependent.variable.name = varname[j], num.trees = 5000,
                      importance = "impurity_corrected")
#        print(mod.imp)
        ## select statistically significant taxa (initial list of candidates)
        imp0 <- importance_pvalues(mod.imp)
        imp[,j] <- imp0[,2]
        namesav <- dimnames(imp0)[[1]][imp0[,2] < 0.05]

        mod <- ranger(data = ss0[, c(varname[j], namesav)],
                      dependent.variable.name = varname[j], num.trees = 5000,
                      importance = "permutation")

        predsav[incvec,j] <- mod$predictions
#        print(mod)

        if (! is.null(keytaxa[[j]][1])) {

            ## some exploratory work to see what's going on in
            ## individual trees
            runrpart <- F
            if (runrpart) {

                ## slashes don't work in rpart, so change them to
                ## underscores
                for (k in 1:length(keytaxa))
                    keytaxa[[k]] <- gsub("/", "_", keytaxa[[k]])

                ## randomly select mtry taxa for running one tree
                ## closer look at the tree structure
                set.seed(13)
                nit <- 100
                require(rpart)
                for (kk in 1:nit) {
                    taxasamp <- c(sample(keytaxa[[j]], 6), "ENCHYTRAEIDAE")
                    formstr <- paste(varname[j], taxasamp[1], sep = "~")
                    for (i in 2:length(taxasamp))
                        formstr <- paste(formstr, taxasamp[i], sep = "+")
                    mod <- rpart(as.formula(formstr), data = ss0, method = "anova",
                                 control = list(minbucket = 5))
                    print(mod$variable.importance)

                    a <- which("ENCHYTRAEIDAE" == names(mod$variable.importance))
                    if (length(a) > 0) {
                        if (a < 5) {
                            dev.new()
                            print(mod)
                            plot(mod)
                            text(mod)
                        }
                    }
                }
                stop()

                require(mgcv)
                resp1 <- ss0$LEPIDOSTOMA == 1 & ss0$PROTOPTILA == 0
                mod1 <- gam(resp1 ~ s(ss0[, varname[j]], k = 4),
                            family = "binomial")
                resp2 <- ss0$LEPIDOSTOMA == 1 & ss0$PROTOPTILA == 1
                mod2 <- gam(resp2 ~ s(ss0[, varname[j]], k = 4),
                            family = "binomial")

                pdf(width = 9, height = 6, pointsize = 10,
                    file = "plots.pdf")
                par(mar = c(4,4,3,1), mfrow = c(2,3), bty = "l")

                iord <- order(ss0[, varname[j]])
                for (i in 1:length(keytaxa[[j]])) {
                    resp3 <- ss0[, keytaxa[[j]][i]] == 1
                    mod3 <- gam(resp3 ~ s(ss0[, varname[j]], k = 4),
                                family = "binomial")
                    predout3 <- predict(mod3, type = "response")

                    plot(ss0[iord, varname[j]], predout3[iord], ylim = c(0,1),
                         type = "l", axes = F)
                    logtick.exp(0.001, 10, c(1),  c(F,F))
                    axis(2)
                    mtext(keytaxa[[j]][i], side = 3, line= 0)
                }
                dev.off()
                stop()
            }

            ## rerun model with just key taxa and
            ## save predictions over previous ones
            modk <- ranger(data = ss0[, c(varname[j], keytaxa[[j]])],
                      dependent.variable.name = varname[j], num.trees = 5000,
                      importance = "permutation")
#            print(modk)
            predsav[incvec,j] <- modk$predictions
            print(summary(lm(modk$predictions ~ ss0[, varname[j]])))
        }


#        print(namesav)

        ## rerun model only with significant taxa
#        mod <- ranger(data = ss0[, c(varname[j], namesav)],
#                      dependent.variable.name = varname[j], num.trees = 2000,
#                      importance = "permutation")
#        print(mod)

        ## plot inferred vs observed env conditions
        predplot <- T
        if (predplot) {
            plot(modk$predictions, ss0[, varname[j]], pch = 21,
                 col = "grey", bg = "white",
                 xlab = paste(lab0[j], "(Predicted)"),
                 ylab = paste(lab0[j], "(Observed)"), axes = F)
            mtext(paste("R2 =", round(modk$r.squared, digits = 2)),
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
        }
        ## otherwise plot important taxa
        else {
            require(pdp)

            ## try different numbers of taxa to optimize model performance
            testtaxa <- F
            if (testtaxa) {
                ## calculate partial dependence relationship
                ## for each taxon selected by importance significance
                cat("Number selected taxa:", length(namesav), "\n")
                for (i in 1:length(namesav)) {
                    cat(i, " ")
                    if (floor(i/10) == i/10) cat("\n")
                    flush.console()
                    pout <- partial(mod, pred.var = namesav[i], plot = FALSE)
                    peff[namesav[i],j] <- (pout$yhat[2] - pout$yhat[1])/
                        diff(range(predsav[,j], na.rm = T))
                }
                cat("\n")

                ## select TRUE here to try different cutoff values
                ## for peff to maximize model performance
                trycuts <- F
                if (trycuts) {
                    cutoffs <- c(0.001,0.003,  0.007,  0.011,  0.015, 0.017,
                                 0.019, 0.021, 0.023)
                    pred.err <- rep(NA, times = length(cutoffs))
                    ntaxa <- rep(NA, times = length(cutoffs))
                    for (k in 1:length(cutoffs)) {
                        selvec <- abs(peff[,j]) > cutoffs[k]
                        selvec[is.na(selvec)] <- F
                        namesnew <- dimnames(peff)[[1]][selvec]
                        ntaxa[k] <- length(namesnew)
                        if (length(namesnew) > 10) {
                            modloc <- ranger(data = ss0[, c(varname[j], namesnew)],
                                             dependent.variable.name = varname[j],
                                             num.trees = 5000,
                                             importance = "permutation")
                            pred.err[k] <- modloc$prediction.error
                        }
                    }
                    print(pred.err)
                    print(ntaxa)
                    dev.new()
                    plot(cutoffs, pred.err)
                    stop()
                }
                else {
                    selvec <- abs(peff[,j]) > cutsav[j]
                    selvec[is.na(selvec)] <- F
                    keytaxa[[j]] <- dimnames(peff)[[1]][selvec]
                }
            }

            ## plot taxa assuming keytaxa is known
            doplot <- T
            if (doplot) {
                peff <- rep(NA, times = length(keytaxa[[j]]))
                names(peff) <- keytaxa[[j]]
                for (i in 1:length(keytaxa[[j]])) {
                    cat(i, " ")
                    if (floor(i/10) == i/10) cat("\n")
                    flush.console()
                    pout <- partial(modk, pred.var = keytaxa[[j]][i], plot = FALSE)
                    peff[i] <- (pout$yhat[2] - pout$yhat[1])/
                        diff(range(predsav[,j], na.rm = T))
                }
                cat("\n")

                incvec <- peff < 0
                peff.neg <- peff[incvec]

                iord <- order(peff.neg)
                plot(peff.neg[iord], 1:length(peff.neg), axes = F,
                     xlim = range(c(0, peff.neg)),
                     xlab = paste("Change in", lab0[j]), ylab = "")
                abline(v = 0, lty = "dashed")
                axis(1)
                axis(2, at = 1:length(peff.neg), lab = names(peff.neg)[iord],
                     cex.axis = 0.65)
                box(bty = "l")

                incvec <- peff > 0
                peff.pos <- peff[incvec]

                iord <- order(peff.pos)
                plot(peff.pos[iord], 1:length(peff.pos), axes = F,
                     xlim = range(c(0, peff.pos)),
                     xlab = paste("Change in", lab0[j]), ylab = "")
                abline(v = 0, lty = "dashed")
                axis(1)
                axis(2, at = 1:length(peff.pos), lab = names(peff.pos)[iord],
                     cex.axis = 0.65)
                box(bty = "l")
                dev.off()
                stop()
            }

#            plot(peff[,3], peff[,4])

        }

    }
    dev.off()

    return(predsav)
#    return(keytaxa)

    dev.new()
    pairs(predsav)
    print(cov(predsav, use = "pair"))
    return()
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

#keytaxa<-explore2(ss, site.data)
predsav.rf <- explore2(ss, site.data, keytaxa)
