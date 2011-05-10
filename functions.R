print_OUT<-function(string,LOG){
  print(paste(date(),"      ",string),quote = FALSE)  
}

lambda<-function(p){ 
  median( qchisq(1-p,df=1),na.rm=T  )/0.456
}
apply_gc<-function(p,lambda){
    pchisq( qchisq(1-p,df=1)/lambda,df=1,lower.tail=F)
}

my_make_positive_definite <- function (cov_in) {
  cov_out<-NULL
  if ( is.positive.definite(cov_in) == FALSE){  
    cov_out<-make.positive.definite(cov_in)
  }
  if (is.null(cov_out) ==FALSE){
    if(is.positive.definite(cov_out) == FALSE){
      cov_out = cov_in;
      diag(cov_out) <- 1.0001; 
    }
    if(is.positive.definite(cov_out) == FALSE){
      diag(cov_out) <- 1.001; 
    }
    if(is.positive.definite(cov_out) == FALSE){
      diag(cov_out) <- 1.01; 
    }
  } else {
    cov_out<-cov_in  
  }
  return(cov_out)
}

modified_fisher<-function(p,w=rep(1/length(p),length(p)),cor = diag(1,length(p))){
	diag(cor) <- 0;
  cor = (3.263*abs(cor) + 0.710*(abs(cor)^2) + 0.027*(abs(cor)^3))
  df = 8/( sum( w * w * cor ) + 4*sum(w^2) )
  chi_stat = (sum(-2 * log(p) * w,na.rm=TRUE)/2) * df;
  return(list('p' = pchisq(chi_stat,df=df,lower.tail=F), "chi"=chi_stat,"df"=df))
}

gates_transform_corr<-function(x){
    x<-0.2982*(x^6) - 0.0127*(x^5) + 0.0588*(x^4) + 0.0099*(x^3) + 0.6281*(x^2) - 0.0009*x; 
    return (x);
}

z_fix_and_random_effects<-function(z,w = rep(1/length(z),length(z)),cov = diag(1,length(z))){
    z_fix<- sum(z*w)/sum(w,na.rm=TRUE)
    v_fix<- sum( w %x% t(w) * cov ,na.rm=TRUE)
    
    Q = 0.0
    {
  	  cor = cov
  		diag(cor) <- 0;
		  cor = (3.263*abs(cor) + 0.710*(abs(cor)^2) + 0.027*(abs(cor)^3))
		  df = 8/(sum( w %x% t(w) * cov ,na.rm=TRUE) + 4*sum(w^2,na.rm=TRUE))
		  Q_naive = sum (w * ((z_fix - z)^2),na.rm=TRUE);
		  Q = qchisq( pchisq(df*0.5*Q_naive, df=df), df= length(z) - 1)
      if (is.na(Q) == T) { Q = 0.0;}
	  }
    I_squared = 0.0;
    if (Q > (length(z) - 1)){
		  I_squared = 100*(Q - (length(z) - 1))/Q
	  }
    tau_squared = 0.0;
    if (Q > (length(z) - 1)){
		  tau_squared = (Q - (length(z) - 1))/(sum(w,na.rm=TRUE) - sum(w^2,na.rm=TRUE)/sum(w,na.rm=TRUE))
	  }
    w_star = 1/( tau_squared + 1/w);
    z_random<- sum(z*w_star,na.rm=TRUE)/sum(w_star,na.rm=TRUE)
    v_random<- sum( w_star %x% t(w_star) * cov ,na.rm=TRUE)
    
    return(
        list( "z_fix" = z_fix,
              "v_fiz" = v_fix,
              "P_FIX" = pnorm(z_fix,sd=sqrt(v_fix),lower.tail=F),
              "z_random" = z_random,
              "v_random" = v_random,
              "P_RANDOM" = pnorm(z_random,sd=sqrt(v_random),lower.tail=F),
              "Q" = Q,
              "I2" = I_squared,
              "tau_squared" = tau_squared
        )
      )
}

rmvnorm<-function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
    method = c("eigen", "svd", "chol")) {
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
        warning("sigma is numerically not symmetric")
    }
    method <- match.arg(method)
    if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
            t(ev$vectors)
    }
    else if (method == "svd") {
        sigsvd <- svd(sigma)
        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    else if (method == "chol") {
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[, o]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
}

report_advance<-function(index,report,tag) {
  if (index > 0){
   if ( ( (index /report ) - round(index/report) ) == 0) {
     print_OUT(paste("   '->Done with [ ",index," ] ",tag,sep="")); 
   }
  }
}