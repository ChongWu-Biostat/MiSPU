# library(CompQuadForm)
MiSPU = function(y, X,tree,cov = NULL, resample = c("perm", "sim", "boot"), model = c("gaussian", "binomial"), pow = c(2:8, Inf), n.perm = 1000){
        
    if (any(is.na(y))){
        ids = which(is.na(y))
        stop(paste("subjects", ids, "has missing response, please remove before proceed \n"))
    }
    
    
    if(is.null(X)==FALSE){
        if(NROW(X)!= length(y)) stop("Dimensions of X and y don't match.")
    }
    
    if (!is.null(X)){
        if (any(is.na(X))){
            stop("NAs in  covariates X, please impute or remove subjects which has missing covariates values")
        }
    }
    
    ## calcualte cum and br.len (branch length)
    GuniF.cum = GUniFrac_cum(X,tree)
    cum = GuniF.cum$cum
    br.len = GuniF.cum$br.len
    br.len = as.matrix(br.len)

    # create generalized taxa proportions
    cum2 = cum *matrix(rep(br.len,each = dim(cum)[2]),nrow = dim(cum)[1],ncol = dim(cum)[2],byrow = TRUE)
    
    # applying SPAT and aSPAT
    aSPU.res <- aSPU(Y, t(cum2), cov = cov, resample = resample,model = model, pow = pow, n.perm = n.perm)

    return(aSPU.res)
}
