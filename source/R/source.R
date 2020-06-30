
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
  if (!is.matrix(data) && !is.data.frame(data)) stop("ERROR: data must be a matrix or data frame")
  data <- as.matrix(data)
  if (!is.matrix(data)) stop("ERROR: data cannot be coerced to a matrix")  
  if (!is.numeric(data)) stop("ERROR: data must be numeric")
  if (any(data < 0)) stop("ERROR: data must be non-negative") 

  data

}

check_op <- function(op, nc, nr) {

  valid <- c("min.rank", "max.rank", "k.fold", "n.seeds",
             "max.iter", "em.eps", "plot", "print",
             "seeds", "kfold.vec", "min.value", "get.summary",
             "n.cores")
  def   <- list(1, 10, 10, 30,
                2000, 1e-5, TRUE, 1,
                NULL, NULL, 0.0001, 1,
                1)
  op  <- default.list(op, valid, def)
  nm  <- names(list)
  tmp <- !(nm %in% valid)
  if (any(tmp)) {
    err <- nm[tmp]
    str <- paste(err, collapse=", ", sep="")
    msg <- paste("ERROR: the option(s) ", str, " are not valid", sep="")
    stop(msg)
  }
  if (op$min.rank < 1) stop("ERROR with option min.rank")
  if (op$max.rank < 1) stop("ERROR with option max.rank")
  if (op$min.rank > op$max.rank) stop("ERROR with option min.rank and/or max.rank")
  if (op$k.fold < 2) stop("ERROR with option k.fold")
  if (op$n.seeds < 1) stop("ERROR with option n.seeds")
  if (op$max.iter < 1) stop("ERROR with option max.iter")
  if ((op$em.eps <= 0) || (op$em.eps > 1)) stop("ERROR with option em.eps")
  kvec <- op[["kfold.vec", exact=TRUE]]
  if (!length(kvec)) {
    op$kfold.vec <- 1:(op$k.fold)
  } else {
    tmp <- !(kvec %in% 1:(op$k.fold))
    if (any(tmp)) stop("ERROR with option kfold.vec") 
  }
  nseeds <- op$n.seeds
  seeds  <- op[["seeds", exact=TRUE]]
  mseeds <- length(seeds)
  if (mseeds && (mseeds != nseeds)) op$n.seeds <- mseeds
  if (!mseeds) op$seeds <- floor(runif(nseeds, min=1, max=1e8))

  m <- min(nr, nc)
  if (op$min.rank > m) stop(paste("ERROR: option min.rank cannot exceed ", m, sep=""))
  op$max.rank <- min(m, op$max.rank)

  op <- set_op_par(op)
  op$algorithm <- 1
  
  op

}

set_op_par <- function(op) {

  n.cores <- op[["n.cores", exact=TRUE]]
  if (!length(n.cores)) n.cores <- 1
  os <- .Platform$OS.type
  if (tolower(os) == "windows") n.cores <- 1
  rvec    <- (op$min.rank):(op$max.rank)
  n.runs  <- length(op$seeds)*length(op$kfold.vec)*length(rvec)

  if ((n.cores < 2) || (n.runs < 2)) {
    op$n.cores <- 1
    #return(op)
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
  for (i in 1:n.cores) {
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

get_idx_mat <- function(NR, NC, k, Kfold) {

  ret <- matrix(data=FALSE, nrow=NR, ncol=NC)
  for(COL in 1:NC) {
    idx <- seq((COL+k-1)%%Kfold, NR, by=Kfold)
    if (idx[1] == 0) idx <- idx[-1]
    ret[idx, COL] <- TRUE
  }
  ret

}

loglike <- function(input, input.k.hat, idxMat, delta_denom) {

  mat1 <- input.k.hat    
  tmp  <- mat1 < delta_denom
  if (any(tmp)) mat1[tmp] <- delta_denom
  vec  <- input[idxMat]
  vec1 <- input.k.hat[idxMat]
  lik  <- vec*(log(vec) - log(mat1[idxMat])) 
  tmp  <- vec == 0
  if (any(tmp)) lik[tmp] <- 0
  lik <- lik - vec + vec1

  ret <- sum(lik)

  ret
}

loglike_new <- function(inputVec, logInputVec, input.k.hatVec, delta_denom) {

  # logInputVec should be finite even if inputVec = 0
  lvec      <- log(input.k.hatVec)
  tmp       <- input.k.hatVec < delta_denom
  lvec[tmp] <- log(delta_denom)
  lik       <- inputVec*(logInputVec - lvec) - inputVec + input.k.hatVec
  ret       <- sum(lik)

  ret
}


get_iargs <- function(x, rank, op=NULL, rvec=NULL) {

  if (is.null(op)) {
    ret <- c(nrow(x), ncol(x), rank, 0, 0, rank, 0, 0)
  } else {
    ret <- c(nrow(x), ncol(x), rank, op$k.fold, 
            op$max.iter, max(rvec), op$print, length(rvec))
  }
  if (length(ret) != 8) stop("ERROR with iargs")
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


my_nmf_C <- function(x, rank, seed) {

  iargs <- get_iargs(x, rank)
  dargs <- get_dargs()
  set.seed(seed)

  nr    <- nrow(x)
  nc    <- ncol(x)
  retW  <- matrix(-1, nrow=nr, ncol=rank)
  retH  <- matrix(-1, nrow=rank, ncol=nc)

  tmp   <- .C("C_call_nmf", as.numeric(x), as.integer(iargs), as.numeric(dargs), 
               retW=as.numeric(retW), retH=as.numeric(retH), PACKAGE="SUITOR")
  retW  <- matrix(tmp$retW, nrow=nr, ncol=rank, byrow=FALSE)
  retH  <- matrix(tmp$retH, nrow=rank, ncol=nc, byrow=FALSE)
  mat   <- retW %*% retH
  
  list(mat=mat, W=retW, H=retH)

}

call_nmf <- function(x, rank, seed, alg=1) {

  if (alg) {
    ret <- my_nmf_C(x, rank, seed)
  } else {
    #tmp <- nmf(x, rank=rank, nrun=1, seed=seed)
    #mat <- tmp@fit@W %*% tmp@fit@H
    #ret <- list(mat=mat, W=tmp@fit@W, H=tmp@fit@H)
  }

  ret

}

EM_alg <- function(input, input.k, idxMat, rank, seed, train.adj, 
                   minValue=0.0001, maxiter=100, EPS=1e-5, alg=1) {

  conv <- 0
  iter <- 0
  while(1) {
    iter <- iter + 1
    ## update NMF
    tmp  <- input.k < minValue
    if (any(tmp)) input.k[tmp] <- minValue
    input.k.hat <- call_nmf(input.k, rank, seed, alg=alg)$mat
    delta_denom <- min(input.k.hat[input.k.hat != 0])/2
    lik1        <- loglike(input, input.k.hat, !idxMat, delta_denom)
    lik1        <- sqrt(lik1/train.adj)
    if ((iter > 1) && (abs((lik0-lik1)/lik0) < EPS)) {
      conv <- 1
      break
    }
    if(iter >= maxiter) break
    lik0            <- lik1 
    input.k[idxMat] <- input.k.hat[idxMat]
  }

  list(converged=conv, niter=iter, loglike=lik1, input.k.hat=input.k.hat,
       delta_denom=delta_denom)

}

ECM_alg <- function(input, input.k, W0, H0, idxMat, rank, seed, train.adj, 
                   minValue=0.0001, maxiter=100, EPS=1e-5, alg=1) {

  conv <- 0
  iter <- 0
  nr   <- nrow(input)
  nc   <- ncol(input)
  mat1 <- matrix(1, nr, nc)

  inputVec         <- input[!idxMat]
  logInputVec      <- log(inputVec)
  tmp              <- inputVec == 0
  logInputVec[tmp] <- 0

  while(1) {
    iter <- iter + 1
    ## update NMF
    tmp  <- input.k < minValue
    if (any(tmp)) input.k[tmp] <- minValue  

    tmp         <- t(W0)
    H1          <- H0*(tmp %*% (input.k/(W0 %*% H0)))/(tmp %*% mat1)
    tmp         <- t(H1)
    W1          <- W0*((input.k/(W0 %*% H1)) %*% tmp)/(mat1 %*% tmp)
    input.k.hat <- W1 %*% H1

    delta_denom <- min(input.k.hat[input.k.hat != 0])/2


    lik1        <- loglike_new(inputVec, logInputVec, input.k.hat[!idxMat], delta_denom)
    #lik1       <- loglike(input, input.k.hat, !idxMat, delta_denom)
    lik1        <- sqrt(lik1/train.adj)

    if ((iter > 1) && (abs((lik0-lik1)/lik0) < EPS)) {
      conv <- 1
      break
    }
    if(iter >= maxiter) break
    lik0            <- lik1 
    W0              <- W1
    H0              <- H1
    input.k[idxMat] <- input.k.hat[idxMat]
  }

  list(converged=conv, niter=iter, loglike=lik1, input.k.hat=input.k.hat,
       delta_denom=delta_denom)

}

suitor_inner <- function(r, seed, input, input.k, idxMat, minValue, maxiter, EPS, print, alg=1) {

   NR              <- nrow(input)
   NC              <- ncol(input)
   tmp             <- call_nmf(input.k, r, seed, alg=alg) 
   input.k.hat     <- tmp$mat
   W               <- tmp$W
   H               <- tmp$H
   tmp             <- NULL
   delta_denom     <- min(input.k.hat[input.k.hat != 0])/2
   train.adj       <- NR*NC - sum(idxMat)   

   #lik0            <- loglike(input, input.k.hat, !idxMat, delta_denom) 
   #lik0            <- sqrt(lik0/train.adj)
   input.k[idxMat] <- input.k.hat[idxMat]

   tmp <- try(ECM_alg(input, input.k, W, H, idxMat, r, seed, train.adj, minValue=minValue, 
                      maxiter=maxiter, EPS=EPS, alg=alg), silent=TRUE)
   if ("try-error" %in% class(tmp)) tmp <- list(converged=FALSE) 
   conv <- tmp$converged
   if (!conv) {
      warning("EM algorithm did not converge")
      err.tmp   <- NA
      train.tmp <- NA 
   } else {
      if (print) cat(paste("EM algorithm converged in ", tmp$niter, " iterations\n", sep=""))
      input.k.hat <- tmp$input.k.hat
      delta_denom <- tmp$delta_denom

      tmp         <- NULL
      # input.k is not used and re-initialized at beginning of loop
      err.tmp      <- loglike(input, input.k.hat, idxMat, delta_denom) 
      train.tmp    <- loglike(input, input.k.hat, !idxMat, delta_denom) 
   }

   list(test=err.tmp, train=train.tmp)

}

suitor_main <- function(input, op) {

  n <- op$n.cores
  if (n < 2) {
    ret <- suitor_seq_C(input, op, op$parMat)
  } else {
    clus <- makeForkCluster(n)
    #clus  <- makeCluster(n)
    registerDoParallel(clus)
    tmp <- suitor_par(input, op)
    stopCluster(clus)

    # Combine results
    nruns <- nrow(op$parMat)
    cx    <- colnames(tmp[[1]])  
    ret   <- matrix(data=NA, nrow=nruns, ncol=length(cx))
    colnames(ret) <- cx
    b <- 0  
    for (i in 1:length(tmp)) {
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

  funs <- c("get_idx_mat", "get_iargs_nmf", "get_dargs_nmf", "loglike",
            "my_nmf_C", "call_nmf", "suitor_inner", "suitor_par_main", 
            "ECM_alg")
  i <- -1
  #foreach(i=1:n, .verbose=FALSE, .inorder=FALSE) %dopar% {
  #  suitor_par_main(input, op, mat[a[i]:b[i], , drop=FALSE])  
  #}
  foreach(i=1:n, .verbose=FALSE, .inorder=FALSE) %dopar% {
    suitor_seq_C(input, op, mat[a[i]:b[i], , drop=FALSE])  
  }

}

initReturnMat <- function(n) {

  tmp           <- c("Rank", "k", "Seed", "Error.Train", "Error.Test")
  ret           <- matrix(data=NA, nrow=n, ncol=length(tmp))  
  colnames(ret) <- tmp 

  ret
}

suitor_par_main <- function(input, op, parMat) {

  Kfold    <- op$k.fold
  minValue <- op$min.value
  maxiter  <- op$max.iter
  EPS      <- op$em.eps
  print    <- op$print
  alg      <- op$algorithm
  n        <- nrow(parMat)
  NR       <- nrow(input)
  NC       <- ncol(input)

  ret      <- initReturnMat(n)
  k0       <- -1
  
  for (i in 1:n) {
     v    <- parMat[i, ]
     r    <- v[1]
     k    <- v[2]
     seed <- v[3]

     if (k != k0) {
       input.k         <- input
       idxMat          <- get_idx_mat(NR, NC, k, Kfold)
       input.k[idxMat] <- NA
 
       for (mut in 1:NR) {
          tmp               <- is.na(input.k[mut,])
          input.k[mut, tmp] <- median(input.k[mut,!tmp])
       }
      
       ## initialize W and H
       tmp             <- input.k < minValue
       input.k[tmp]    <- minValue
     }
     k0       <- k
     tmp      <- suitor_inner(r, seed, input, input.k, idxMat, minValue, maxiter, EPS, print, alg=1)
     ret[i, ] <- c(r, k, seed, tmp$train, tmp$test)
  }

  ret

}

suitor_seq_C <- function(input, op, parMat) {

  seeds  <- unique(parMat[, 3])
  nseeds <- length(seeds)
  N      <- nrow(parMat)
  ret    <- initReturnMat(N)
  MISS   <- -9999.0e100
  MISS2  <- -9999.0e99
 
  b <- 0
  for (i in 1:nseeds) {
    seed  <- seeds[i]
    set.seed(seed) 
    tmp   <- parMat[, 3] == seed
    rvec  <- parMat[tmp, 1]
    kvec  <- parMat[tmp, 2]
    nk    <- length(kvec)
    iargs <- get_iargs(input, 0, op=op, rvec=rvec)
    dargs <- get_dargs(op=op)
    err1  <- as.numeric(rep(MISS, nk))
    err2  <- as.numeric(rep(MISS, nk)) 
   
    tmp  <- .C("C_call_suitor", as.numeric(input), as.integer(iargs), 
               as.numeric(dargs), as.integer(rvec), as.integer(kvec), 
                ret_train=err1, ret_test=err2, PACKAGE="SUITOR")
    err1 <- tmp$ret_train
    err2 <- tmp$ret_test
    tmp  <- err1 < MISS2
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
  }  

  ret

}

suitor_seq <- function(input, op) {

  rankVec  <- (op$min.rank):(op$max.rank)
  KfoldVec <- op$kfold.vec
  Kfold    <- op$k.fold
  minValue <- op$min.value
  maxiter  <- op$max.iter
  EPS      <- op$em.eps
  seeds    <- op$seeds
  print    <- op$print
  alg      <- op$algorithm
  NR       <- nrow(input)
  NC       <- ncol(input)
   
  nseeds   <- length(seeds)
  nranks   <- length(rankVec)
  ret      <- initReturnMat(nseeds*nranks*length(KfoldVec))

  row <- 0
  for(r in rankVec){
    if (print) cat(paste("rank = ", r, "\n", sep=""))
    for(k in KfoldVec){
      if (print > 1) cat(paste("k = ", k, "\n", sep=""))

      input.k <- input

      ############################
      ## divide training/test elements in V
      idxMat          <- get_idx_mat(NR, NC, k, Kfold)
      input.k[idxMat] <- NA
 
      ############################
      ## Initial guess for test elements
      for (mut in 1:NR) {
        tmp               <- is.na(input.k[mut,])
        input.k[mut, tmp] <- median(input.k[mut,!tmp])
      }
 
 
   
      tmp             <- input.k < minValue
      input.k[tmp]    <- minValue
      input.k0        <- input.k

      for (i in 1:nseeds) {
        seed            <- seeds[i]
        if (print > 2) cat(paste("seed = ", seed, "\n", sep=""))


        input.k         <- input.k0
        tmp             <- call_nmf(input.k, r, seed, alg=alg) 
        input.k.hat     <- tmp$mat

        W               <- tmp$W
        H               <- tmp$H
        tmp             <- NULL
        delta_denom     <- min(input.k.hat[input.k.hat != 0])/2
        train.adj       <- NR*NC - sum(idxMat)   

        #lik0            <- loglike(input, input.k.hat, !idxMat, delta_denom) 
        #lik0            <- sqrt(lik0/train.adj)
        input.k[idxMat] <- input.k.hat[idxMat]

        tmp <- try(ECM_alg(input, input.k, W, H, idxMat, r, seed, train.adj, minValue=minValue, 
                      maxiter=maxiter, EPS=EPS, alg=alg), silent=TRUE)
        if ("try-error" %in% class(tmp)) tmp <- list(converged=FALSE) 
        conv <- tmp$converged
        if (!conv) {
          warning("EM algorithm did not converge")
          err.tmp   <- NA
          train.tmp <- NA 
        } else {
          if (print > 2) cat(paste("EM algorithm converged in ", tmp$niter, " iterations\n", sep=""))
          input.k.hat <- tmp$input.k.hat
          delta_denom <- tmp$delta_denom

          tmp         <- NULL
          # input.k is not used and re-initialized at beginning of loop

          err.tmp      <- loglike(input, input.k.hat, idxMat, delta_denom) 
          train.tmp    <- loglike(input, input.k.hat, !idxMat, delta_denom) 
        }
        row            <- row + 1
        ret[row, ]     <- c(r, k, seed, train.tmp, err.tmp)
      }
    }
  }
  
  list(all.results=ret)

}

getSummary <- function(obj, NC, NR=96) {

  opt.r  <- NA
  M      <- NR*NC
  ranks  <- obj[, 1]
  kvec   <- obj[, 2]
  uranks <- sort(unique(ranks))
  nranks <- length(uranks)
  uk     <- sort(unique(kvec))
  nk     <- length(uk)
  train  <- obj[, 4]
  test   <- obj[, 5]
  CV.tr  <- rep(NA, nranks)
  CV.te  <- rep(NA, nranks)

  cx     <- c("Rank", "Type", "MSErr", paste("fold", uk, sep=""))
  tab    <- matrix(data="", nrow=2*nranks, ncol=length(cx))
  colnames(tab) <- cx

  tmp0 <- is.finite(train) & is.finite(test)
  row  <- 0
  for (i in 1:nranks) {
    rank <- uranks[i]
    tmpr <- ranks == rank
    CV.1 <- rep(NA, nk)
    CV.2 <- CV.1
    for (j in 1:nk) {
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
    mse1       <- sqrt(CV.tr[i]/M)
    mse2       <- sqrt(CV.te[i]/M)
    row        <- row + 1
    tab[row, ] <- c(rank, "Train", mse1, CV.1)
    row        <- row + 1
    tab[row, ] <- c(rank, "Test", mse2, CV.2)
  }

  tab <- as.data.frame(tab, stringsAsFactors=FALSE)
  tmp <- (1:ncol(tab))[-2]
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

  ranks <- as.numeric(x[, "Rank"])
  cls   <- c("blue", "green")
  min.r <- min(ranks, na.rm=TRUE)
  max.r <- max(ranks, na.rm=TRUE)
  mse   <- as.numeric(x[, "MSErr"])
  min.e <- min(mse, na.rm=TRUE)
  max.e <- max(mse, na.rm=TRUE)
  plot(0, 0, type="n", xlab="Rank", ylab="MSE", xlim=c(min.r, max.r),
       ylim=c(min.e, max.e))
  tmp <- x[, "Type"] %in% "Train"
  points(ranks[tmp], mse[tmp], type="l", col=cls[1])
  tmp <- x[, "Type"] %in% "Test"
  v1  <- ranks[tmp]
  v2  <- mse[tmp]
  points(v1, v2, type="l", col=cls[2])
  j   <- which.min(v2)
  r   <- v1[j]
  e   <- v2[j]
  points(r, e, type="p", pch=19, col="red")
  leg <- c("Train", "Test")
  legend("top", leg, fill=cls, horiz=TRUE)

  NULL
}

# Function to assign a default value to an element in a list
default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} # END: default.list

getSeqsFromList <- function(inlist) {

  ncomb <- 1
  nc    <- length(inlist)
  for (i in 1:nc) ncomb <- ncomb*length(inlist[[i]])
  nn    <- names(inlist)

  mat <- matrix(NA, nrow=ncomb, ncol=nc)
  if (length(nn)) colnames(mat) <- nn
  for (j in 1:nc) {
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

####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
std.my.update.h <- function(v, w, h, nbterms=0L, ncterms=0L, copy=TRUE)
{	
	.Call("my_update_H", v, w, h, nbterms, ncterms, copy)
}
R_std.my.update.h <- function(v, w, h, wh=NULL)
{	
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# divergence-reducing NMF iterations
	# H_au = H_au ( sum_i [ W_ia V_iu / (WH)_iu ] ) / ( sum_k W_ka ) -> each row of H is divided by a the corresponding colSum of W
	h * crossprod(w, v / wh) / colSums(w)	
}
std.my.update.w <- function(v, w, h, nbterms=0L, ncterms=0L, copy=TRUE)
{	
	.Call("my_update_W", v, w, h, nbterms, ncterms, copy)
}
R_std.my.update.w <- function(v, w, h, wh=NULL)
{			
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# W_ia = W_ia ( sum_u [H_au A_iu / (WH)_iu ] ) / ( sum_v H_av ) -> each column of W is divided by a the corresponding rowSum of H
	#x2 <- matrix(rep(rowSums(h), nrow(w)), ncol=ncol(w), byrow=TRUE); 
	#w * tcrossprod(v / wh, h) / x2;
	sweep(w * tcrossprod(v / wh, h), 2L, rowSums(h), "/", check.margin = FALSE) # optimize version?
	
}

my_update.brunet_R <- function(i, v, w, h, eps=.Machine$double.eps)
{
	
	# standard divergence-reducing NMF update for H
	h <- R_std.my.update.h(v, w, h)

	# standard divergence-reducing NMF update for W
	w <- R_std.my.update.w(v, w, h)
	
	#every 10 iterations: adjust small values to avoid underflow 
	if( i %% 10 == 0 ){
		h[h<eps] <- eps
		w[w<eps] <- eps
	}
	
	#return the modified model	
	return(list(W=w, H=h))
	
}

my_update.brunet <- function(i, v, w, h, eps=.Machine$double.eps)
{
	
	# standard divergence-reducing NMF update for H	
	h <- std.my.update.h(v, w, h, nbterms=0L, ncterms=0L)
	
	# standard divergence-reducing NMF update for W
	w <- std.my.update.w(v, w, h, nbterms=0L, ncterms=0L)

	
	#every 10 iterations: adjust small values to avoid underflow
	if( i %% 10 == 0 ){
		h[h<eps] <- eps
		w[w<eps] <- eps
	}
	
	list(W=w, H=h)
	
}

myStopFun <- function(consold, inc, h, stopconv=40) {
				
   # construct connectivity matrix
   index <- apply(h, 2, function(x) which.max(x) )  # max in each col
   cons  <- outer(index, index, function(x,y) ifelse(x==y, 1,0))

   changes <- cons != consold
   if( !any(changes) ) {
     inc <- inc + 1 # connectivity matrix has not changed: increment the count
   } else {
      consold <- cons
      inc     <- 0   # else restart counting
   }
									
   # assume convergence is connectivity stops changing 
   if (inc > stopconv) {
     ret <- TRUE
   } else {
     ret <- FALSE
   }

  list(rc=ret, consold=consold, inc=inc)
		
}


my_nmf_R <- function(x, rank, seed, maxiter=2000, check.interval=10, 
                   eps.update=.Machine$double.eps) {

  set.seed(seed)
  nr      <- nrow(x)
  nc      <- ncol(x)
  m       <- max(x, na.rm=TRUE)
  W       <- matrix(runif(nr*rank, min=0, max=m), nrow=nr,   ncol=rank, byrow=FALSE)
  H       <- matrix(runif(nc*rank, min=0, max=m), nrow=rank, ncol=nc,   byrow=FALSE)
  p       <- ncol(x)
  consold <- matrix(0, p, p)
  inc     <- 0
  conv    <- 0
  i       <- 0
  x0      <- x


  while(1) {		
    i <- i + 1
			
    # update the matrices
    tmp     <- my_update.brunet(i, x, W, H, eps=eps.update)
    W       <- tmp$W
    H       <- tmp$H
    tmp     <- NULL

    # test convergence only every 10 iterations
    if (i %% check.interval == 0) {
      tmp     <- myStopFun(consold, inc, H, stopconv=40) 
      conv    <- tmp$rc
      consold <- tmp$consold
      inc     <- tmp$inc
    }

    # if the strategy ask for stopping, then stop the iteration
    if (conv || i >= maxiter ) break	
    			
  }

  W %*% H

}
