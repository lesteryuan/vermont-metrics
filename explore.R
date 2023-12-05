## 12.4.2023
## initial exploration

explore <- function() {

    bug.data$bio_sample <- factor(bug.data$bio_sample)

    ## replace text N/A with NA
    incvec <- bug.data$Genus == "N/A"
    bug.data$Genus[incvec] <- NA
    incvec <- bug.data$GenusGroup == "N/A"
    bug.data$GenusGroup[incvec] <- NA
    incvec <- bug.data$SubFamilyOrTribe == "N/A"
    bug.data$SubFamilyOrTribe[incvec] <- NA

    taxon <- bug.data$Genus
    print(sum(is.na(taxon)))

    incvec <- is.na(bug.data$Genus)
    taxon[incvec] <- bug.data$GenusGroup[incvec]
    print(sum(is.na(taxon)))
    incvec <- is.na(taxon)
    taxon[incvec] <- bug.data$SubFamilyOrTribe[incvec]
    print(sum(is.na(taxon)))

    print(nrow(bug.data))
    incvec <- is.na(taxon)
    print(sum(incvec))

    numocc <- table(taxon)
    numocc.sav <- numocc[numocc > 30]
    tnames <- names(numocc.sav)

    ## drop slash groups
    w <- regexpr("/", tnames)
    tnames <- tnames[w == -1]

    ## drop uid
    tnames <- tnames[tnames != "UID" & tnames != "UNID" & tnames != "GROUP"]

    ## drop genus
    w <- regexpr("GENUS", toupper(tnames))
    tnames <- tnames[w == -1]

    ## drop thie group
    w <- regexpr("GROUP", tnames)
    tnames <- tnames[w==-1]

    ## make site-species matrix
    ss <-matrix(FALSE, ncol = length(tnames),
                nrow = length(levels(bug.data$bio_sample)))
    dimnames(ss) <- list(levels(bug.data$bio_sample), tnames)
    for (i in 1:length(tnames)) {
        incvec <- taxon == tnames[i]
        incvec[is.na(incvec)] <- F
#        z <- tapply(bug.data$density_m2[incvec], bug.data$bio_sample[incvec],
#                    sum)
                                        #        z[is.na(z)] <- 0
              #  ss[names(z), tnames[i]] <- as.vector(z)
        y <- unique(bug.data$bio_sample[incvec])
        ss[y, tnames[i]] <- TRUE
        ss[,tnames[i]]<-factor(ss[,tnames[i]])
    }

    ## spot check SS
#    sitep <- levels(bug.data$bio_sample)[10]
#    incvec <- bug.data$bio_sample == sitep
#    print(bug.data[incvec,c ("Genus", "density_m2")])
#    print(ss[sitep,ss[sitep,]>0])

    print(dim(ss))
    ss <- data.frame(ss)
    ss$bio_sample <- levels(bug.data$bio_sample)

    site.data$bio_sample <- factor(site.data$bio_sample)
    ## drop 1 duplicate bio_sample
    site.data <- site.data[unique(site.data$bio_sample),]

    ss <- merge(ss, site.data,by = "bio_sample")
    print(dim(ss))

    ss$conductivity <- log(ss$conductivity)
    ss$tp_ug <- log(ss$tp_ug)
    ss$nitrate_mg <- log(ss$nitrate_mg)
    incvec <- ss$turbidity == 0
    incvec[is.na(incvec)] <- F
    minval <- 0.5*min(ss$turbidity[!incvec], na.rm = T)
    ss$turbidity[incvec] <- minval
    ss$turbidity <- log(ss$turbidity)

    png(width = 6, height =  10, units = "in", res = 600,
        file = "taxop.png")
#    par(mar = c(4,4,3,2), mfrow = c(1,2), mgp = c(2.3,1,0))
   par(mar = c(4,8,3,2),las=1, mfrow = c(1,2), mgp = c(2.3,1,0))
    require(ranger)
    varname <- c("turbidity", "silt_rating")
    lab0 <- c("Turbidity", "Silt rating")
    logt <- c(T, F)
    for (j in 1:length(varname)) {
        incvec <- ! is.na(ss[, varname[j]])
        print(sum(incvec))
        ss0 <- ss[incvec,]

        mod <- ranger(data = ss0[, c(varname[j], tnames)],
                      dependent.variable.name = varname[j], num.trees = 2000,
                      importance = "impurity_corrected")
        print(mod)
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
        imp <- importance_pvalues(mod)
        namesav <- dimnames(imp)[[1]][imp[,2] < 0.01]
        print(namesav)

        mod <- ranger(data = ss0[, c(varname[j], namesav)],
                      dependent.variable.name = varname[j], num.trees = 2000,
                      importance = "permutation")
        print(mod)
        peff <- rep(NA, times = length(namesav))
        names(peff) <- namesav

        require(pdp)
        for (i in 1:length(namesav)) {
            pout <- partial(mod, pred.var = namesav[i], plot = FALSE)
            peff[i] <- pout$yhat[2] - pout$yhat[1]
        }
        print(peff)
        iord <- order(peff)
        plot(peff[iord], 1:length(peff), axes = F,
             xlab = paste("Change in", lab0[j]), ylab = "")
        abline(v = 0, lty = "dashed")
        axis(1)
        axis(2, at = 1:length(peff), lab = namesav[iord], cex.axis = 0.65)
        box(bty = "l")

    }
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


explore()
#cluster.env(site.data)
