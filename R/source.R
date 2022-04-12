
suitor <- function(data, op=NULL) {

    data <- check_data(data)
    nr   <- nrow(data)
    nc   <- ncol(data)
    op   <- check_op(op, nc, nr)
    ret  <- suitor_main(data, op)
    if (op$get.summary) {
        ret <- getSummary(ret$all.results, nc, NR=nr)
        if (op$plot) plotErrors(ret$summary) 
    }
    ret$op <- op
    ret
}

check_data <- function(data) {

    if (!length(data)) stop("ERROR with input data")
    if (!is.matrix(data) && !is.data.frame(data)) 
        stop("ERROR: data must be a matrix or data frame")
    data <- as.matrix(data)
    if (!is.matrix(data)) stop("ERROR: data cannot be coerced to a matrix")  
    if (!is.numeric(data)) stop("ERROR: data must be numeric")
    if (any(data < 0)) stop("ERROR: data must be non-negative") 
    if (any(!is.finite(data))) stop("ERROR: data must be non-negative") 

    data

}

check_op_valid <- function(op) {

    valid <- c("min.rank", "max.rank", "k.fold", "n.starts",
        "max.iter", "em.eps", "plot", "print",
        "kfold.vec", "min.value", "get.summary",
        "n.cores", "type")
    def   <- list(1, 10, 10, 30,
        2000, 1e-5, TRUE, 1,
        NULL, 0.0001, 1,
        1, NULL)
    op  <- default.list(op, valid, def)
    nm  <- names(op)
    tmp <- !(nm %in% valid)
    if (any(tmp)) {
        err <- nm[tmp]
        str <- paste(err, collapse=", ", sep="")
        msg <- paste("ERROR: the option(s) ", str, " are not valid", sep="")
        stop(msg)
    }
    op

}

check_op <- function(op, nc, nr, which=1) {

    # which  1=suitor, 2=extractWH

    op <- check_op_valid(op)

    if (which == 2) {
        op$k.fold    <- 1
        op$kfold.vec <- 1
    }
    if (op$min.rank < 1) stop("ERROR with option min.rank")
    if (op$max.rank < 1) stop("ERROR with option max.rank")
    if (op$min.rank > op$max.rank) 
        stop("ERROR with option min.rank and/or max.rank")
    if ((which == 1) && (op$k.fold < 2)) stop("ERROR with option k.fold")
    if (op$n.starts < 1) stop("ERROR with option n.starts")
    if (op$max.iter < 1) stop("ERROR with option max.iter")
    if ((op$em.eps <= 0) || (op$em.eps > 1)) stop("ERROR with option em.eps")
        kvec <- op[["kfold.vec", exact=TRUE]]
    if (!length(kvec)) {
        op$kfold.vec <- seq_len(op$k.fold)
    } else {
        tmp <- !(kvec %in% seq_len(op$k.fold))
        if (any(tmp)) stop("ERROR with option kfold.vec") 
    }

    op$seeds <- seq_len(op$n.starts)

    m   <- min(nr, nc)
    msg <- paste("ERROR: option min.rank cannot exceed ", m, sep="")
    if (op$min.rank > m) stop(msg)
    op$max.rank <- min(m, op$max.rank)

    op <- set_op_par(op)
    op$algorithm <- 1

    op$type <- set_op_type(op[["type", exact=TRUE]])
    type    <- op$type
    if ((length(type) != 1) || !is.character(type)) 
        stop("ERROR with option type")

    op

}

set_op_type <- function(type) {

    if (is.null(type)) {
        os <- tolower(.Platform$OS.type)
        if ("unix" %in% os) {
            type <- "FORK"
        } else {
            type <- "PSOCK" 
        }
    }
    type <- toupper(removeWhiteSpace(type))
    type

}

set_op_par <- function(op) {

    n.cores <- op[["n.cores", exact=TRUE]]
    if (!length(n.cores)) n.cores <- 1
    rvec    <- (op$min.rank):(op$max.rank)
    n.runs  <- length(op$seeds)*length(op$kfold.vec)*length(rvec)

    if ((n.cores < 2) || (n.runs < 2)) {
        op$n.cores <- 1
    }

    mat       <- getSeqsFromList(list(op$seeds, op$kfold.vec, rvec))
    tmp       <- mat[, 1]
    mat[, 1]  <- mat[, 3]
    mat[, 3]  <- tmp
    op$parMat <- mat
    if (nrow(mat) != n.runs) stop("ERROR")
    m         <- floor(n.runs/n.cores)
    rem       <- n.runs - m*n.cores
    parStart  <- rep(NA, n.cores)
    parEnd    <- rep(NA, n.cores)
    b         <- 0
    for (i in seq_len(n.cores)) {
        a           <- b + 1
        b           <- a + m - 1
        if (rem) {
            b   <- b + 1
            rem <- rem - 1
        }
        parStart[i] <- a
        parEnd[i]   <- b
    }
    if (parEnd[n.cores] != n.runs) stop("ERROR2")
    op$parStart <- parStart
    op$parEnd   <- parEnd
    op$n.cores  <- n.cores

    op

}

get_iargs <- function(x, rank, op=NULL, rvec=NULL, startNum=0) {

    if (is.null(op)) {
        ret <- c(nrow(x), ncol(x), rank, 0, 0, rank, 0, 0, startNum)
    } else {
        ret <- c(nrow(x), ncol(x), rank, op$k.fold, 
            op$max.iter, max(rvec), op$print, length(rvec), startNum)
    }
    if (length(ret) != 9) stop("ERROR with iargs")
    ret
}

get_dargs <- function(op=NULL) {

    if (is.null(op)) {
        ret <- c(.Machine$double.eps, 0, 0)
    } else {
        ret <- c(.Machine$double.eps, op$em.eps, op$min.value)
    }
    if (length(ret) != 3) stop("ERROR with dargs")
    ret
}

my_nmf_C <- function(x, rank, op=NULL) {

    iargs <- get_iargs(x, rank)
    dargs <- get_dargs(op=op)
    nr    <- nrow(x)
    nc    <- ncol(x)
    retW  <- matrix(-1, nrow=nr, ncol=rank)
    retH  <- matrix(-1, nrow=rank, ncol=nc)

    if (length(op)) {
        minv <- op$min.value
        x    <- as.numeric(x)
        tmp  <- x < minv
        tmp[is.na(tmp)] <- FALSE
        if (any(tmp)) x[tmp] <- minv 
    }

    tmp   <- .C("C_call_nmf", as.numeric(x), as.integer(iargs), 
        as.numeric(dargs), retW=as.numeric(retW), 
        retH=as.numeric(retH), PACKAGE="SUITOR")
    retW  <- matrix(tmp$retW, nrow=nr, ncol=rank, byrow=FALSE)
    retH  <- matrix(tmp$retH, nrow=rank, ncol=nc, byrow=FALSE)
    mat   <- retW %*% retH

    list(mat=mat, W=retW, H=retH)

}

suitor_main <- function(input, op) {

    n <- op$n.cores
    if (n < 2) {
        ret <- suitor_seq_C(input, op, op$parMat)
    } else {
        if (op$type == "FORK") {
            clus  <- makeForkCluster(n)
        } else {
            clus <- makeCluster(n, type=op$type)
            obj  <- c("suitor_seq_C", "initReturnMat", "get_iargs", "get_dargs")
            clusterExport(cl=clus, obj, envir=environment())
        }   
        registerDoParallel(clus)
        tmp <- suitor_par(input, op)
        stopCluster(clus)

        # Combine results
        nruns <- nrow(op$parMat)
        cx    <- colnames(tmp[[1]])  
        ret   <- matrix(data=NA, nrow=nruns, ncol=length(cx))
        colnames(ret) <- cx
        b <- 0  
        for (i in seq_len(length(tmp))) {
            mat <- tmp[[i]]
            a   <- b + 1
            b   <- a + nrow(mat) - 1
            if (b > nruns) stop("ERROR combining results")
            ret[a:b, ] <- mat
        }
    }

    list(all.results=ret)

}

suitor_par <- function(input, op) {

    mat  <- op$parMat
    a    <- op$parStart
    b    <- op$parEnd
    n    <- op$n.cores
    i    <- -1

    if (op$type == "FORK") {
        ret <- foreach(i=seq_len(n), .verbose=FALSE, .inorder=FALSE) %dopar% {
            suitor_seq_C(input, op, mat[a[i]:b[i], , drop=FALSE])  
        }
    } else {
        ret <- foreach(i=seq_len(n), .verbose=FALSE, .inorder=FALSE, 
            .packages='SUITOR') %dopar% {
            suitor_seq_C(input, op, mat[a[i]:b[i], , drop=FALSE])  
        }
    }
    ret
}

initReturnMat <- function(n) {

    tmp           <- c("Rank", "k", "Start", "Error.Train", "Error.Test",
        "EM.niter")
    ret           <- matrix(data=NA, nrow=n, ncol=length(tmp))  
    colnames(ret) <- tmp 

    ret
}

suitor_seq_C <- function(input, op, parMat) {

    seeds       <- unique(parMat[, 3])
    nseeds      <- length(seeds)
    N           <- nrow(parMat)
    ret         <- initReturnMat(N)
    MISS        <- -9999.0e100
    MISS2       <- -9999.0e99

    b <- 0
    for (i in seq_len(nseeds)) {
        seed     <- seeds[i]
        tmp      <- parMat[, 3] == seed
        rvec     <- parMat[tmp, 1]
        kvec     <- parMat[tmp, 2]
        nk       <- length(kvec)
        iargs    <- get_iargs(input, 0, op=op, rvec=rvec, startNum=i)
        dargs    <- get_dargs(op=op)
        err1     <- as.numeric(rep(MISS, nk))
        err2     <- as.numeric(rep(MISS, nk)) 
        convVec  <- as.integer(rep(0, nk))
        niterVec <- as.integer(rep(0, nk))

        tmp  <- .C("C_call_suitor", as.numeric(input), as.integer(iargs), 
            as.numeric(dargs), as.integer(rvec), as.integer(kvec), 
            ret_train=err1, ret_test=err2, ret_conv=convVec, 
            ret_niter=niterVec, PACKAGE="SUITOR")
        err1     <- tmp$ret_train
        err2     <- tmp$ret_test
        niterVec <- tmp$ret_niter
        tmp      <- err1 < MISS2
        if (any(tmp)) err1[tmp] <- NA
        tmp  <- err2 < MISS2
        if (any(tmp)) err2[tmp] <- NA

        a           <- b + 1
        b           <- a + nk - 1
        tmp         <- a:b
        ret[tmp, 1] <- rvec
        ret[tmp, 2] <- kvec
        ret[tmp, 3] <- seed
        ret[tmp, 4] <- err1
        ret[tmp, 5] <- err2
        ret[tmp, 6] <- niterVec
    }  

    ret

}

getSummary <- function(obj, NC, NR=96) {

    M      <- NR*NC
    ranks  <- obj[, 1]
    kvec   <- obj[, 2]
    uranks <- sort(unique(ranks))
    nranks <- length(uranks)
    uk     <- sort(unique(kvec))
    Kfold  <- max(uk, na.rm=TRUE)
    nk     <- length(uk)
    train  <- obj[, 4]
    test   <- obj[, 5]
    CV.tr  <- CV.te  <- rep(NA, nranks)
    tab    <- matrix(data="", nrow=2*nranks, ncol=3+nk)
    tmp0   <- is.finite(train) & is.finite(test)
    row    <- 0
    for (i in seq_len(nranks)) {
        rank <- uranks[i]
        tmpr <- ranks == rank
        CV.1 <- rep(NA, nk)
        CV.2 <- CV.1
        for (j in seq_len(nk)) {
            k    <- uk[j]
            tmpk <- kvec == k
            tmp  <- tmpr & tmpk & tmp0
            tmp[is.na(tmp)] <- FALSE
            vec1 <- train[tmp]
            vec2 <- test[tmp]
            if (!length(vec1)) next
            arg  <- which.min(vec1)
            if (length(arg)) {
                CV.1[j] <- vec1[arg]
                CV.2[j] <- vec2[arg]
            }
        }
        CV.tr[i]   <- sum(CV.1, na.rm=TRUE)
        CV.te[i]   <- sum(CV.2, na.rm=TRUE)
        mse1       <- sqrt(CV.tr[i]/(M*(Kfold-1)))
        mse2       <- sqrt(CV.te[i]/M)
        row        <- row + 1
        tab[row, ] <- c(rank, "Train", mse1, CV.1)
        row        <- row + 1
        tab[row, ] <- c(rank, "Test", mse2, CV.2)
    }

    return(getSummary_setRetObj(tab, CV.te, uranks, obj, uk))

}

getSummary_setRetObj <- function(tab, CV.te, uranks, obj, uk) {

    cx     <- c("Rank", "Type", "MSErr", paste("fold", uk, sep=""))
    colnames(tab) <- cx
    opt.r <- NA
    tab   <- as.data.frame(tab, stringsAsFactors=FALSE)
    tmp   <- (seq_len(ncol(tab)))[-2]
    for (i in tmp) tab[, i] <- as.numeric(tab[, i])

    tmp  <- is.finite(CV.te)
    vec1 <- CV.te[tmp]
    vec2 <- uranks[tmp]
    if (length(vec1)) {
        j     <- which.min(vec1)
        opt.r <- vec2[j]
    }

    list(rank=opt.r, summary=tab, all.results=obj)

}

plotErrors <- function(x) {

    errors <- as.numeric(x$MSErr)
    ranks  <- as.numeric(x$Rank)
    types  <- x$Type

    tmp     <- types %in% "Test"
    idx.min <- which.min(errors[tmp])
    x.p     <- ranks[tmp][idx.min]
    y.p     <- errors[tmp][idx.min]

    tmp <- ggplot(data=x, aes(x=ranks, y=errors, group=types)) +
        geom_line(aes(color=types),size=1)+
        scale_x_continuous(breaks=seq(1, 15, 1),name="Number of signatures")+
        scale_y_continuous(name="Prediction error")+
        scale_colour_manual(name  ="",
            breaks=c("Test", "Train"),
            labels=c("Validation", "Training"),
            values=c("red", "blue"))+
        geom_point(aes(x=x.p, y=y.p), colour="red")+  
        theme_bw()+
        theme(legend.position="bottom")

    print(tmp)

    NULL

}

# Function to assign a default value to an element in a list
default.list <- function(inList, names, default) {

    # inList      List
    # names       Vector of names of items in inList
    # default     List of default values to assign if a name is not found
    #             The order of default must be the same as in names.

    n1 <- length(names)
    n2 <- length(default)
    if (n1 != n2) stop("ERROR: in calling default.list")

    if (is.null(inList)) inList <- list()

    listNames <- names(inList)
    for (i in seq_len(n1)) {
        if (!(names[i] %in% listNames)) {
            inList[[names[i]]] <- default[[i]]   
        }
    }

    inList

} # END: default.list

getSeqsFromList <- function(inlist) {

    ncomb <- 1
    nc    <- length(inlist)
    for (i in seq_len(nc)) ncomb <- ncomb*length(inlist[[i]])
    nn    <- names(inlist)

    mat <- matrix(NA, nrow=ncomb, ncol=nc)
    if (length(nn)) colnames(mat) <- nn
    for (j in seq_len(nc)) {
        vec  <- inlist[[j]]
        nvec <- length(vec)
        if (j == 1) {
            m <- ncomb/nvec
        } else {
            m <- m/nvec
        }
        mat[, j] <- rep(vec, each=m)
    }

    mat

}

removeWhiteSpace <- function(str, leading=1, trailing=1) {

    if ((leading) && (trailing)) {
        ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
    } else if (leading) {
        ret <- gsub("^\\s+", "", str, perl=TRUE)
    } else if (trailing) {
        ret <- gsub("\\s+$", "", str, perl=TRUE)
    } else {
        ret <- str
    }

    ret

} # END: removeWhiteSpace

