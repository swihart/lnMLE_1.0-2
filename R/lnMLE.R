#.packageName <- "lnMLE"
#.First.lib <- function(lib, pkg) library.dynam("lnMLE", pkg, lib)


logit.normal.mle <- function(meanmodel, logSigma, id, n=NULL, beta=NULL, alpha=NULL, model="marginal", lambda=0.0, r = 20, maxits=50, tol = 1e-3, data = sys.frame(sys.parent()) )
{
#
# maximum likelihood for logistic-normal model
# random intercept std dev design
#
	mint <- match( model, c("marginal","conditional") )
	if( is.na(mint) ){
	  stop("specify model as 'marginal' or 'conditional'.")
	}

	mrerror <- match( r, c(3,5,10,20,50) )
	if( is.na(mrerror) ){
	  stop("Error: please choose r = 3, 5, 10, 20, or 50")
	}


# Get x terms and y:	
  mf <- model.frame(meanmodel, data)
  terms <- attr(mf, "terms")
  tl.x <- attr(terms, "term.labels")
  intr <- attr(terms, "intercept")
  if(intr > 0) tl.x <- c("(Intercept)", tl.x)
  y <- model.response(mf, "numeric")
  nobs <- length(y)
  x <- model.matrix(meanmodel, mf, contrasts)
  	
  if( is.null(n) ) n<-rep( 1, length(y) )
  nclust<-length( table( id ) )

# Get z terms:	
  mf <- model.frame(logSigma, data)
  terms <- attr(mf, "terms")
  tl.z <- attr(terms, "term.labels")
  intr.z <- attr(terms, "intercept")
  if(intr.z > 0) tl.z <- c("(Intercept)", tl.z)
  z <- model.matrix(logSigma, mf, contrasts)
	
  if(is.null(beta)) {
    # Provide a reasonable estimate for beta:
	fitbeta <- glm(y ~ -1 + x, family=binomial)
	beta <- fitbeta$coef
	
  } else {
  
    if(length(beta) != length(x[1,])) {stop("Number of estimates in beta is ", length(beta), 
			" should equal ", length(x[1,]), " (number ofcovariates in x).")
	}
  }
  p <- length(x[1,])
  bnames <- tl.x
  if( is.null( bnames ) ) bnames <- paste("beta",1:p)
  
  q <- length(z[1,])
  anames <- tl.z
  if( is.null( anames ) ) anames <- paste("alpha",1:q)
 
  if(is.null(alpha)) {
  
	# Provide a reasonable initial value for alpha:
	alpha <- rep(0.1, q)
	
  } else {
    if(length(alpha) != length(z[1,])) {stop("Number of estimates in alpha is ", length(alpha), 
			" should equal ", length(z[1,]), " (number ofcovariates in z).")
    }
  }
  
  		
	mod.cov<-matrix( 0, p+q, p+q )

	flag<-0
#
	#options(contrasts = c("contr.treatment", "contr.poly" ))
	beta.m <- glm( cbind( y, n-y ) ~ -1 + x, family=binomial )$coef
  	logL <- 0.0
	int.parms <- c( p, q, r, nobs, nclust )
	eta <- as.vector( x%*%beta.m )
#
	z <- .C("logit_normal_mle",
		as.double(id),
		as.double(y),
		as.double(n),
		as.double(x),
		beta = as.double(beta),
		as.double(z),
		alpha = as.double(alpha),
		as.integer(mint),
		eta = as.double(eta),
		as.integer(int.parms),
		as.double(lambda),
		logL = as.double(logL),
		iter=as.integer(maxits),
		tol=as.double(tol),
		mod.cov=as.double(mod.cov),
		as.double(beta.m),
		flag=as.integer(flag),
		PACKAGE = "lnMLE" )	
#
	cov.model =  matrix(z[["mod.cov"]],p+q,p+q)
	beta<-as.vector(z[["beta"]])
	names(beta)<-bnames
	pbeta<-length(beta)
	se.beta<-sqrt(diag(cov.model))[1:pbeta]
	betas<-cbind( beta, se.beta, beta/se.beta )
	dimnames( betas )<-list( names(beta), c("estimate","std. err.","  Z  ") )

	alpha<-as.vector(z[["alpha"]])
	names(alpha)<-anames
	qalpha<-length(alpha)
	se.alpha<-sqrt(diag(cov.model))[(pbeta+1):(pbeta+qalpha)]
	alphas<-cbind( alpha, se.alpha, alpha/se.alpha )
	dimnames( alphas )<-list( names(alpha), c("estimate","std. err.","  Z  ") )
#
out <- list( 	title = "ML Estimation for Logistic-Normal Models",
		dispersion = "logistic-normal random effects",
		model = model,
		beta = beta,
		betas = betas,
		alpha = alpha,
		alphas = alphas, 
		eta = z[["eta"]],
		cov.model =  cov.model,
		logL = z[["logL"]],
		options=list( lambda=lambda, r=r ),
		tol=z[["tol"]], 
		iter=z[["iter"]],
		flag=z[["flag"]] )
#
class(out) <- "logit.normal.mle"
out
}

"print.logit.normal.mle" <-
function( x, digits=NULL,... ){
#
if (is.null(digits)) 
        digits <- options()$digits
else options(digits = digits)

cat("\n")
cat( paste(" ",x$title ) )
cat("\n")
cat( paste("\n  model =", x$model ) )
cat("\n")
cat(paste("\n  options:     lambda =", round(x$options$lambda,3) ) )
cat(paste("\n                    r =", x$options$r, "\n"))
#beta<-fit$beta
#p<-length(beta)
#se.beta<-sqrt(diag(x$cov))[1:p]
#betas<-cbind( beta, se.beta, beta/se.beta )
#dimnames( betas )<-list( names(x$beta), c("estimate","std. err.","  Z  ") )
cat("\n\n  Mean Parameters: \n\n")
print( x$betas, digits=digits)
#alpha<-fit$alpha
#q<-length(alpha)
#se.alpha<-sqrt(diag(fit$cov))[(p+1):(p+q)]
#tab<-cbind( alpha, se.alpha, alpha/se.alpha )
#dimnames( tab )<-list( names(fit$alpha), c("estimate","std. err.","  Z  ") )
cat("\n\n  Variance Components: \n\n")
print( x$alphas, digits=digits)
cat(paste("\n\n  Maximized logL =", round(x$logL, digits=digits),"\n\n"))
invisible(x)
#end...
}

