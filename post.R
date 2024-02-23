## post process pairwise importance values
post <- function() {
    load("peff.single.rda")
    peff.sing <- peff
    print(peff.sing)

    ## additive effect
    peff.add <- matrix(NA, nrow = length(peff.sing), ncol = length(peff.sing))
    dimnames(peff.add)[[1]] <- names(peff.sing)
    dimnames(peff.add)[[2]] <- names(peff.sing)
    for (i in 2:length(peff.sing)) {
        for (j in 1:(i-1)) {
            peff.add[i,j] <- abs(peff.sing[i]) + abs(peff.sing[j])
        }
    }
    print(summary(peff.add))

    load("peff.rda")
    peff[is.na(peff)] <- 0
    peff[is.infinite(peff)] <- 0

    peff.sc <- peff/peff.add
    incvec <- peff.sc == 0
    incvec[is.na(incvec)] <- F
    peff.sc[incvec] <- NA
    print(summary(peff.sc))
    hist(peff.sc)
    print(summary(as.vector(peff.sc)))

    r0 <- rank(as.vector(peff.sc), na.last = "keep")
    print(length(r0))
    dim(r0) <- c(nrow(peff), ncol(peff))
    dimnames(r0) <- dimnames(peff)
    print(r0)

    ntot <- max(r0, na.rm = T)
#    for (i in ntot:(ntot-10)) {
    for (i in 1:10) {
        selvec <- r0 == i
        selvec[is.na(selvec)] <- F
        dim(selvec) <- c(nrow(peff), ncol(peff))
        colnum <-which( apply(selvec, 2, any, na.rm = T))
        rownum <- which( apply(selvec, 1, any, na.rm = T))

        cat(dimnames(r0)[[2]][colnum], dimnames(r0)[[1]][rownum], sep = "/")
        cat("\n")
        cat(peff.sc[rownum, colnum])
        cat("\n")

        print(peff[rownum, colnum])
        print(peff.sing[rownum])
        print(peff.sing[colnum])

    }

}

post()
