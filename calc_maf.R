calc_maf <- function(x)
{
    nr <- nrow(x)     
    nc <- ncol(x)     
    
    x[(x!=0) & (x!=1) & (x!=2)] <- NA
    x <- as.matrix(x)
    
    ## calc_n
    n0 <- apply(x == 0, 1, sum,na.rm=T)
    n1 <- apply(x == 1, 1, sum,na.rm=T)
    n2 <- apply(x == 2, 1, sum,na.rm=T)
    
    n <- n0 + n1 + n2
    
    ## calculate allele frequencies
    p <- ((2*n0)+n1)/(2*n)
    q <- 1 - p
    maf <- pmin(p, q)

    # MODIFIED 21 Oct 2012:  prior to this version, we had "mono=(mgf<0)" instead of "mono<(maf<0)"
    res <- data.frame( n=n, n0=n0, n1=n1, n2=n2, p=p, maf=maf,
                       mono=(maf<=0), loh=(n1<=0), 
                       stringsAsFactors=F )
    row.names(res) <- row.names(x)
    res
}
