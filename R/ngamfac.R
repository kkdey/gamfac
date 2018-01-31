#' @title Normal Gamma Lasso factor analysis for Gaussian data
#'
#' @description
#' Fit a sparse normal factor analysis model equipped with  Gamma Lasso shrinkage
#' priors on Gaussian noise data.
#'
#' @param X a matrix (features by samples) assumed to comprise of Gaussian distributed random variables.
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
ngamfac <- function (X, numfac, cl = NULL,
                     covars = NULL, numfac_init = NULL, 
                     override=FALSE, NUM_ITER = 5){
  
  ############  checking the X matrix   ####################
  
  library(gamlr)
  library(parallel)

  ###########   checking the covariates    ###################
  if(!is.null(covars)){
    if(nrow(covars) != ncol(X)){
      stop("The first dimension of covariates must be equal to the number of columns of X matrix.")
    }
  }
  m <- nrow(X)
  n <- ncol(X)
  
  NA_IND <- is.na(X)
  mean_X <- rowMeans(X, na.rm=TRUE)
  norm_X <- X - mean_X
  norm_X[NA_IND] <- 0
  
  adjust <- 8

  norm_z <- centerscale(norm_X)
  
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
    argl$family <- "normal"
    argl$x <- v[,-ncol(v)]
    C <- dim(X)[1]
    chunks <- lapply(1:C, 
                     function(i) X[i,])
    
    # cat("distributed run.\n") 
    #  print(cl)
    mods <- parLapply(cl,chunks,onerun_normal,argl=argl) 
    h <- do.call(rbind, mods)
    
    
    D <- dim(X)[2]
    chunks2 <- lapply(1:D, 
                      function(i) X[,i])
    argl2 <- argl
    argl2$x <- h[,1:numfac]
    mods <- parLapply(cl,chunks2,onerun_normal,argl=argl2) 
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

onerun_normal <- function(xj, argl){
  argl$y <- xj
  fit <- do.call(cv.gamlr,argl)
  out <- fit$gamlr$beta[,fit$seg.min]
  return(out)
}

