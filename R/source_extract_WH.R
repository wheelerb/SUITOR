
suitorExtractWH <- function(data, rank, op=NULL){

    data <- check_data(data)
    nr   <- nrow(data)
    nc   <- ncol(data)
    if (!is.list(op)) op <- list()
    op$min.rank <- rank
    op$max.rank <- rank
    op          <- check_op(op, nc, nr, which=2)
    ERR_MAT     <- extractWH_main(data, op)

    # Extracting W and H using the seed number
    ii             <- which.min(ERR_MAT[, 2])
    NMF_C          <- my_nmf_C(as.matrix(data), rank=rank, op=op) 
    W_C            <- NMF_C$W
    D_C            <- colSums(W_C)
    if (length(D_C) > 1) {
        D_C          <- diag(D_C)
    } else {
        dim(D_C)     <- c(1, 1) 
    }
    Wn_C           <- apply(W_C, 2, function(x) x/sum(x))
    colnames(Wn_C) <- paste0("denovo ", LETTERS[seq_len(rank)]) 
    rownames(Wn_C) <- rownames(data)
    H_C            <- NMF_C$H
    Hn_C           <- D_C %*% as.matrix(H_C)
    rownames(Hn_C) <- paste0("denovo ", LETTERS[seq_len(rank)]) 

    return(list(W=Wn_C, H=Hn_C))

}

extractWH_main <- function(input, op) {

    n <- op$n.cores
    if (n < 2) {
        ret <- extractWH_seq_C(input, op, op$parMat)
    } else {
        if (op$type == "FORK") {
            clus  <- makeForkCluster(n)
        } else {
            clus  <- makeCluster(n, type=op$type)
            clusterExport(cl=clus, 
                    c("extractWH_seq_C", "my_nmf_C", "get_iargs", "get_dargs"), 
                    envir=environment())
        }
        registerDoParallel(clus)
        tmp <- extractWH_par(input, op)
        stopCluster(clus)

        # Combine results
        ret <- NULL
        for (i in seq_len(length(tmp))) ret <- rbind(ret, tmp[[i]])     
    }
    colnames(ret) <- c("Start", "Error")

    ret

}

extractWH_seq_C <- function(data, op, parMat) {

    seeds   <- unique(parMat[, 3])
    nseeds  <- length(seeds)
    ret     <- rep(NA, nrow(parMat))
    rank    <- op$min.rank
    nc      <- ncol(data)
    print   <- (op$print > 0) && (op$n.cores == 1)
    logdata <- log(data)
    tmp0    <- data == 0
    tmp0[is.na(tmp0)] <- FALSE
    flag0   <- any(tmp0)

    for (i in seq_len(nseeds)) {
        if (print) {
            cat("\r", paste0("start ", i))
            flush.console() 
        }
        seed        <- seeds[i]
        tmp         <- my_nmf_C(as.matrix(data), rank=rank, op=op)
        W1          <- tmp$W
        H1          <- tmp$H
        input.hat   <- W1 %*% H1
        delta_denom <- min(input.hat[input.hat != 0])/2    

        input.hat2        <- input.hat
        tmp               <- input.hat2 < delta_denom
        tmp[is.na(tmp)]   <- FALSE
        if (any(tmp)) input.hat2[tmp] <- delta_denom
        tmp <- data*(logdata - log(input.hat2))
        if (flag0) tmp[tmp0] <- 0
        tmp    <- tmp - data + input.hat
        ret[i] <- sum(tmp, na.rm=TRUE)
    }

    if (print) cat("\n")
    cbind(seeds, ret)

}

extractWH_par <- function(input, op) {

    mat  <- op$parMat
    a    <- op$parStart
    b    <- op$parEnd
    n    <- op$n.cores
    i    <- -1
    if (op$type == "FORK") {
        ret <- foreach(i=seq_len(n), .verbose=FALSE, .inorder=FALSE) %dopar% {
            extractWH_seq_C(input, op, mat[a[i]:b[i], , drop=FALSE])  
        }
    } else {
        ret <- foreach(i=seq_len(n), .verbose=FALSE, .inorder=FALSE,
            .packages="SUITOR") %dopar% {
            extractWH_seq_C(input, op, mat[a[i]:b[i], , drop=FALSE])  
        }
    }
    ret
}
