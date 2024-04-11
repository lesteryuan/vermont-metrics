## examine within-site variability of vt measurements
## 4.11.2023

withinsite <- function(env.all) {
    require(lme4)
    varp <- c("alkalinity", "turbidity")

    cclist <- as.list(rep(NA, times = 2))
    names(cclist) <- varp
    for (j in varp) {
        incvec <- ! is.na(env.all[, j])
        df1 <- env.all[incvec,]

        nvis <- table(df1$site_id)
        idp <- names(nvis[nvis > 1])

        incvec <- rep(F, times = nrow(df1))
        for (i in idp) incvec <- incvec | df1$site_id == i
        df1 <- df1[incvec,]
        
        incvec <- df1[, j] > 0
        df1 <- df1[incvec,]
    
        mod <- lmer(log(df1[, j]) ~ 1 + (1|site_id), data = df1)

        cclist[[j]] <- coef(mod)$site_id
    }

    dftemp <- data.frame(site_id = row.names(cclist[[1]]),
                         alk = cclist[[1]][,1])
    dftemp <- merge(dftemp,
                    data.frame(site_id = row.names(cclist[[2]]),
                               turb = cclist[[2]][,1]),
                    by = "site_id")
    print(dftemp)


    plot(dftemp$turb, dftemp$alk)
    return(dftemp)
}

env.mn <- withinsite(env.all)
    
