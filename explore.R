## 12.4.2023
## initial exploration

## edit bcnt file and generate otu file
explore <- function() {

    ## loaded CSV files into R as bug.data and site.data.

    ## make sure sample id is a factor
    bug.data$bio_sample <- factor(bug.data$bio_sample)

    ## replace text "N/A" with R NA
    incvec <- bug.data$Genus == "N/A"
    bug.data$Genus[incvec] <- NA
    incvec <- bug.data$GenusGroup == "N/A"
    bug.data$GenusGroup[incvec] <- NA
    incvec <- bug.data$SubFamilyOrTribe == "N/A"
    bug.data$SubFamilyOrTribe[incvec] <- NA
    incvec <- bug.data$Family == "N/A"
    bug.data$Family <- NA

    ## change "sp" in species field to NA
    incvec <- tolower(bug.data$Species) == "sp" |
        bug.data$Species == "sp a" | bug.data$Species == "sp b" |
            bug.data$Species == "" | bug.data$Species == "w/o setae" |
                bug.data$Species == "N/A" | bug.data$Species == "spa" |
                    bug.data$Species == "uid" |
                        bug.data$Species == "imm" |
                            bug.data$Species == "uid wsetae"
    bug.data$Species[incvec] <- NA

    incvec <- bug.data$SpeciesGroup == "N/A" | bug.data$SpeciesGroup == "uid" |
        tolower(bug.data$SpeciesGroup) == "group"
    print(table(bug.data$SpeciesGroup))
    bug.data$SpeciesGroup[incvec] <- NA

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
    bug.data <- bug.data[incvec,]

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
explore2 <- function(ss) {

    tnames <- names(ss)[-1]
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
    png(width = 5, height = 5, units = "in", res = 600, pointsize = 8,
        file = "taxop.png")
#    par(mar = c(4,4,3,1), mfrow = c(2,2), mgp = c(2.3,1,0))
    par(mar = c(4,10,1,1), mfrow = c(1,2), mgp = c(2.3,1,0), las = 1)

    require(ranger)
    varname <- c("turbidity", "silt_rating", "alkalinity", "chloride")
    lab0 <- c("Turbidity", "Silt rating", "Alkalinity", "Chloride")
    logt <- c(T, F, T, T)

    predsav <- matrix(NA, ncol = 4, nrow = nrow(ss))
    for (j in 2) {
        incvec <- ! is.na(ss[, varname[j]])
        print(sum(incvec))
        ss0 <- ss[incvec,]

        ## fit initial RF to predict varname (turbidity and silt)
        mod <- ranger(data = ss0[, c(varname[j], tnames)],
                      dependent.variable.name = varname[j], num.trees = 5000,
                      importance = "impurity_corrected")
        ## select statistically significant taxa
        imp <- importance_pvalues(mod)
        namesav <- dimnames(imp)[[1]][imp[,2] < 0.01]

        print(namesav)

        ## rerun model only with significant taxa
        mod <- ranger(data = ss0[, c(varname[j], namesav)],
                      dependent.variable.name = varname[j], num.trees = 2000,
                      importance = "permutation")
        print(mod)

        predsav[incvec,j] <- mod$predictions

        predplot <- F
        if (predplot) {
            plot(mod$predictions, ss0[, varname[j]], pch = 21,
                 col = "grey", bg = "white",
                 xlab = paste(lab0[j], "(Predicted)"),
                 ylab = paste(lab0[j], "(Observed)"), axes = F)
            mtext(paste("R2 =", round(mod$r.squared, digits = 2)),
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

        plottaxa <- T
        if (plottaxa) {
            peff <- rep(NA, times = length(namesav))
            names(peff) <- namesav

            require(pdp)
            for (i in 1:length(namesav)) {
                pout <- partial(mod, pred.var = namesav[i], plot = FALSE)
                peff[i] <- pout$yhat[2] - pout$yhat[1]
            }
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
        }
    }

#    dev.new()
#    plot(predsav[,3], predsav[,4])
    dev.off()
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


#bcnt.otu <- explore()
#cluster.env(site.data)

explore2(ss)
