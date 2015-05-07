core.set <- function(data.set, preset = NULL, coverage = 99, delta = 0.001,
                     coverage_filename = "Coverage.csv",
                     Temp_file = "Temp.csv",
                     Coreset = "Coreset.csv"){
    if (!is.data.frame(data.set)){
        data.set <- as.data.frame(data.set)
    }
    z <- Sys.time()
    print(z)
    counts <- data.set
    data.set <- data.set
    coverage <- as.numeric(coverage)
    nc <- ncol(counts)
    nr <- nrow(counts)
    delta <- delta
    coverage1 <- NULL
    coverage.table <- NULL
    rnames <- rownames(data.set)
    cnames <- colnames(data.set)

    # the number of variables.
    var.num <- vector()
    ide.num <- rep(0, nr)
    var.num <- apply(counts, 1, function(x){length(unique(x[!is.na(x)]))})
    cat(mean(ide.num/var.num*100), "\n")

    mpe<-sum(var.num - 1)
	
    cat("Max possible entries", mpe, "\n")
    result <- NULL
    result.idx <- NULL
    overlap.score <- function(x){
        tmp1 <- x
        tmp1.table <- table(tmp1)
        idx <- match(tmp1, names(tmp1.table))
        dummy <- tmp1.table[idx]
        na.idx <- which(is.na(idx))
        if (length(na.idx) > 0){
            dummy[na.idx] <- 0
        }
            dummy
    }
    prenum <- 0
    if (!is.null(preset)){
        result <- preset
        prenum <- length(result)
        result.idx <- which(match(colnames(data.set), result)>0)
        coreset <- data.set[, result.idx]
        for (i in 1:ncol(coreset)){
            for (k in 1:ncol(counts)){
                idx1 <- which(coreset[, i] == counts[, k])
                counts[idx1, k] <- NA
            }
        }
        counts <- data.set[,-result.idx]
        ide.num <- apply(coreset, 1, function(x){length(unique(x[!is.na(x)]))})
        coverage1 <- mean(ide.num/var.num*100)
        coverage.table <- rbind(coverage.table, c(0, "preset", coverage1, coverage1))
    }

    for (idx in 1:mpe){
        cat(paste0(idx, "-th iteration starts at"), paste(Sys.time()),"\n")
       	cat("iteration",idx, "\n")
        Sys.time()
	
        not.na.counts <- apply(counts, 2, function(x){
            length(which(!is.na(x)))
        })
        rnames = rownames(counts)
        cnames = colnames(counts)
        candidate <- which(not.na.counts == max(not.na.counts))
        step0 <- NULL
        step0 <- data.frame(apply(counts, 1, overlap.score))
        rownames(step0) <- cnames
        if (length(candidate) == 1){
            step01 <- step0[candidate,]
            overlap <- mean(step01[!is.na(step01)])
            names(overlap) <- cnames[candidate]
        } else {
            rownames(step0) <- cnames
            step01 <- step0[candidate,]
            overlap <- apply(step01, 1, function(x){mean(x[!is.na(x)])})
        }

        select <- which(overlap == min(overlap))
        cat("Selected", length(select), "candidate genes for", idx, "\n")

        if (length(select)==1){
            final.select <- select
            final.select.idx <- which(names(select) == colnames(data.set))
            result <- c(result, names(final.select))
            result.idx <- c(result.idx, final.select.idx)
            rm.idx <- candidate[final.select]
        } else {
            minsel <- list()
            minnum <- vector()
            for (i in 1:length(select)){
                minsel[[i]] <- which(!is.na(counts[,candidate[select[i]]]))
                minnum[i] <- min(var.num[minsel[[i]]])
            }
            minselidx <- which(minnum == min(minnum))
            minlen <- length(minselidx)
            if (minlen == 1){
                final.select <- select[minselidx]
            } else {
                final.select <- sample(select[which(minnum == min(minnum))], 1)
            }
            final.select.idx <- which(names(final.select) == colnames(data.set))
            rm.idx <- candidate[final.select]
            result <- c(result, names(final.select))
            result.idx <- c(result.idx, final.select.idx)
        }
		
        cat("Final selected sample is", names(final.select),"\n")
        coreset <- data.frame(data.set[, result.idx])
        colnames(coreset) <- result
        ide.num <- apply(coreset, 1, function(x){length(unique(x[!is.na(x)]))})
        coverage1 = mean(ide.num/var.num*100)

        cat("Coverage is", coverage1, "%", "\n")
        if (idx == 1 & prenum ==0){
            dy <- coverage1
        } else if (idx == 1 & prenum != 0){
            dy <- coverage1 - as.numeric(coverage.table[1,3])
        } else if ( idx > 1 & prenum != 0){
            dy <- coverage1 - as.numeric(coverage.table[idx, 3])
        } else {
            dy <- coverage1 - as.numeric(coverage.table[(idx-1), 3])
        }
        cat("Difference is ", dy*100, "%", "\n")
        coverage.table = rbind(coverage.table, c(idx, names(final.select), coverage1, dy))
        colnames(coverage.table) = c("Iteration","Sample_name", "Coverage", "Difference")

        if(coverage1 >= coverage | dy < delta){
            break
        } else {
            counts <- data.frame(counts[,-rm.idx])
            nc <- ncol(counts)
            if (prenum == 0){
                for (i in 1:nc){
                    idx1 <- which(coreset[, idx] == counts[, i])
                    counts[idx1, i] = NA
                }
            } else {
                for (i in 1:nc){
                    idx1 <- which(coreset[, prenum + idx - 1] == counts[,i])
                    counts[idx1,i] = NA
                }
            }
        }
        write.csv(coverage.table, file = Temp_file, quote = FALSE)
    }
    coverage.table <- as.data.frame(coverage.table)
    cat("GenoCore selects ",length(result), "element for core sets", "\n")
    cat("Running time is", print(Sys.time() - z), "\n")
    write.csv(coverage.table, file = coverage_filename, quote = FALSE)
    write.csv(coreset, file = Coreset, quote = FALSE)
}
