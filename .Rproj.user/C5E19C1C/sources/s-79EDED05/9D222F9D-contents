#' @title Logistic Gamma Lasso factor analysis for Binomial data
#'
#' @description
#' Fit a sparse logistic factor analysis model equipped with  Gamma Lasso shrinkage
#' priors on binomial data. 
#'
#' @note The matrix \code{X} is a binomial matrix. If it is a matrix of presence - absence 
#' (0/1), the \code{max_n} would be 1, if it is a genotype matrix, 0, 1 or 2 coded, 
#' then \code{max_n} would be 2.  We can also have a general X matrix recording the 
#' number of successes at each cell, with \code{max_n} being a matrix of same dimension as X
#' and recording the total numbe rof trials at each cell. 
#'
#' @param X an integer matrix (presence-absence/genotype/ generic binomial success matrix)
#' @param max_n The integer/matrix recording the number of trials from which matrix X is created.
#' Defaults to an integer corresponding to the maximum integer value of X.
#' @param numfac The number of logistic factors (does not include the intercept).
#' @param cl The cluster process created before running the code.
#' @param covars A matrix of known covariates which is modeled together with with the unknown factros. 
#' Defaults to NULL.  
#' @param numfac_init Initial number of factors in the pre-processing steps before gamma lasso analysis.
#' Chosen to be bigger than \code{\numfac}. 
#' @param override optional boolean to bypass Lanczos bidiagonalization
#' SVD. Usually not advised unless encountering a bug in the SVD code.
#' @param NUM_ITER The number of iterative steps of gamma-lasso estimation of factors. 
#' 
#' @return matrix of logistic factors (F) and corresponding loadings (L). 
#' @export
#' @importFrom corpcor fast.svd
#' @importFrom gamlr cv.gamlr
#' @useDynLib lgamfac
#' 
lgamfac <- function (X, max_n = NULL, numfac, cl = NULL,
                     covars = NULL, numfac_init = NULL, 
                     override=FALSE, NUM_ITER = 5){
  
  ############  checking the X matrix   ####################
  
  library(gamlr)
  library(parallel)
  
  if(is.null(max_n)){max_n = max(X)}
  stopifnot( all(X == floor(X)) ) # OK
  if(length(max_n) == 1){
    max_n_out <- matrix(max_n, dim(X)[1], dim(X)[2])
  }else{
    if(nrow(max_n) != nrow(X) | ncol(max_n) != ncol(X)){
      stop("when max_n is a matrix, it must have same dimensions as X")
    }
  }
  Y = max_n_out - X
  if(min(Y) < 0){stop("max_n must be greater than all elements of X matrix")}
  stopifnot( all(Y == floor(Y))) 
  
  ###########   checking the covariates    ###################
  if(!is.null(covars)){
    if(nrow(covars) != ncol(X)){
      stop("The first dimension of covariates must be equal to the number of columns of X matrix.")
    }
  }
  m <- nrow(X)
  n <- ncol(X)

  
  #############  The initial number of factors (> numfac) ############
  
  if(is.null(numfac_init)){
    numfac_init <- numfac + 5
  }else{
    if(numfac_init < (numfac + 5)) numfac_init <- numfac + 5
  }
  
  ############   mean centering the X matrix    ###############
  
  NA_IND <- is.na(X)
  mean_X <- rowMeans(X, na.rm=TRUE)
  norm_X <- X - mean_X
  norm_X[NA_IND] <- 0
  
  adjust <- 8
  
  #######   SVD (using lfa() functionality)  #################
  
  mysvd <- trunc.svd(norm_X, d=numfac_init, adjust=adjust, tol=1e-13, override=override)
  
  rm(norm_X)
  D <- diag(mysvd$d, numfac_init, numfac_init)
  U <- mysvd$u
  V <- mysvd$v
  rm(mysvd)
  
  # form projection from the first SVD  ##############
  
  z <- U %*% D %*% t(V)
  z <- z + mean_X
  z <- z/max_n_out
  
  ################  projected scores into the [0,1]  interval  ###########
  zmin <- min(z)
  zmax <- max(z)
  z1 <- t(apply(z, 1, function(x) return(2/n + (1 - 4/n)*((x - zmin)/(zmax - zmin)))))
  
  z1 <- log(z1/(1-z1))
  
  norm_z <- centerscale(z1)
  
  ########## Second SVD (from lfa() package)  ######################
  
  v <- trunc.svd(norm_z, d=numfac, adjust=adjust, tol=1e-13, override=override)$v
  
  ###  adding the intercept back
  
  if(!is.null(covars)){
    v <- cbind(v,1)
  }else{
    v <- cbind(v, covars, 1)
  }
  
  #############  Fitting sparse Gamma Lasso model iteratively on the SVD vectors  #########
  
  if(is.null(cl)){
    cl <- makeCluster(10,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
    print(cl)
  }
  
  for(iter in 1:NUM_ITER){
    verb <- 0
    argl <- list()
    argl$family <- "binomial"
    argl$x <- v[,-ncol(v)]
    C <- dim(X)[1]
    chunks <- lapply(1:C, 
                     function(i) X[i,])
    
   # cat("distributed run.\n") 
  #  print(cl)
    mods <- parLapply(cl,chunks,onerun_poisson,argl=argl) 
    h <- do.call(rbind, mods)
    
    
    D <- dim(X)[2]
    chunks2 <- lapply(1:D, 
                      function(i) X[,i])
    argl2 <- argl
    argl2$x <- h[,1:numfac]
    mods <- parLapply(cl,chunks2,onerun_poisson,argl=argl2) 
    vnew <- do.call(rbind, mods)
    if(!is.null(covars)){
      v <- cbind(vnew[,1:numfac],covars, 1)
    }else{
      v <- cbind(vnew[,1:numfac], 1)
    }

    cat("we are at iteration,", iter, "\n")
  }
  
  stopCluster(cl)
  
  if(!is.null(covars)){
    colnames(v) <- c(paste0("fac-", 1:numfac), paste0("covar-", 1:dim(covars)[2]), "intercept")
  }else{
    colnames(v) <- c(paste0("fac-", 1:numfac), "intercept")
  }
  if(!is.null(colnames(X))){rownames(v) <- colnames(X)}
  if(!is.null(rownames(X))){rownames(h) <- rownames(X)}
  
  ll <- list("L" = h,
             "F" = v)
  return(ll)
  
}

onerun_poisson <- function(xj, argl){
  argl$y <- xj
  fit <- do.call(cv.gamlr,argl)
  out <- fit$gamlr$beta[,fit$seg.min]
  return(out)
}

