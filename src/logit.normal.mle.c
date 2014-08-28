/* ---------------------------------------------------------------- */
/*                                                                  */
/*  Maximum Likelihood Estimation                                   */
/*   Logistic-Normal Model                                          */
/*    for Clustered Binary (Binomial) Data                          */
/*                                                                  */
/*                                              P. Heagerty         */
/*                                              97/08/07            */
/*                                                                  */
/* ---------------------------------------------------------------- */

#include "chanmatstruct.h"
#include "chanmatfuns.h"


#define cfree free

#define HALFSTEP 1
#define PI 3.1415927
#define RIDGE 1e-1
#define NSAFE 15
#define ETA_TOLERANCE 1e-5
#define ETA_MAX_ITER  50
#define ETA_EPS 1e-7


#define MARGINAL 1
#define CONDITIONAL 2

extern int split();
extern double  *get_GH_z(), *get_GH_w();
extern double  matmaxabs(), antilogit(), logit();
extern void    deconvolve_GH();

extern MATRIX  *matantilogit(), *matxdiagasvec(), *luinv();

void logit_normal_mle( S_id, S_y, S_n, S_x, S_beta, S_z, S_alpha, S_model,
	      S_eta, S_ints, S_lambda, S_logL,
	      S_maxiter, S_tol, S_modcov, S_beta_0, S_flag )

  double  *S_id, *S_y, *S_n, *S_x, *S_beta, *S_z, *S_alpha, *S_eta;
  double  *S_lambda, *S_tol, *S_modcov, *S_beta_0, *S_logL;
  integer *S_ints, *S_maxiter, *S_flag, *S_model;
{
  MATRIX  **Y, **N, **X, **Z, **Eta, **Mu;
  MATRIX  *idin, *yin, *nin, *xin, *zin, *etain;
  MATRIX  *beta, *alpha, *gamma;
  MATRIX  *mu, *eta, *H, *Hi, *U, *Ui;
  MATRIX  *Di, *E;
  MATRIX  *dgda, *delta;
  MATRIX  *beta_0, *mu_0, *lp;
  double  *z, *w;
  double  dmax, tolerance, lambda;
  double  y_i, n_i;
  double  eta_i, gi, dgdai, ps, pm, p_q_s;
  double   log_like, logPi, Pi_s, logPi_s, logPi_indep;
  double  mu_i, detadb, detads, d11, d12, d22;
  double  dEta_dTheta_j, dEta_dTheta_k, d2Eta_dTheta_jk, step;
  int     c, i, j, k, nobs, nclust, p, q, r, one, count, iter, s, Ni;
  int     converge, maxiter, flag, vlink, model;
  FILE    *wfile;

  one = 1;

  vlink = 2;

  model = *S_model;

  p = *(S_ints+0);
  q = *(S_ints+1);
  r = *(S_ints+2);
  nobs = *(S_ints+3);
  nclust = *(S_ints+4);

  tolerance = *S_tol;
  lambda = *S_lambda;
  maxiter = *S_maxiter;

  wfile = fopen("QEE2.out","w");
  fprintf( wfile, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n" );
  fprintf( wfile, " p = %i, q = %i, nobs = %i, nclust = %i\n", p, q, nobs, 
	  nclust );
  fprintf( wfile, "              lambda = %f\n", lambda );
  fprintf( wfile, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n" );
  fclose( wfile );


  from_S( S_beta, &p, &one, beta );
  make_permanent( beta );

  from_S( S_alpha, &q, &one, alpha );
  make_permanent( alpha );

  from_S( S_beta_0, &p, &one, beta_0 );
  make_permanent( beta_0 );
 
  from_S( S_id, &nobs, &one, idin );  
  from_S( S_y, &nobs, &one, yin );
  from_S( S_n, &nobs, &one, nin );
  from_S( S_x, &nobs, &p, xin );
  from_S( S_z, &nobs, &q, zin );
  from_S( S_eta, &nobs, &one, etain );

  Y = (MATRIX **)calloc( nclust, (unsigned)sizeof( struct matrix) );
  N = (MATRIX **)calloc( nclust, (unsigned)sizeof( struct matrix) );
  X = (MATRIX **)calloc( nclust, (unsigned)sizeof( struct matrix) );
  Z = (MATRIX **)calloc( nclust, (unsigned)sizeof( struct matrix) );
  Eta = (MATRIX **)calloc( nclust, (unsigned)sizeof( struct matrix) );
  Mu = (MATRIX **)calloc( nclust, (unsigned)sizeof( struct matrix) );

  make_permanent( idin );

  split( yin, idin, Y );
  split( nin, idin, N );
  split( xin, idin, X );
  split( zin, idin, Z );
  split( etain, idin, Eta );

  destroy_matrix( idin );

  for( c=0; c<nclust; c++ ){
    make_permanent( Y[c] );
    make_permanent( N[c] );
    make_permanent( X[c] );
    make_permanent( Z[c] );
    make_permanent( Eta[c] );
    Mu[c] = create_matrix( Y[c]->nrows, 1, EPHEMERAL );
  }


  /* get quadrature values */

  z = get_GH_z( r );

  w = get_GH_w( r );

  iter=0;
  converge=0;

wfile=fopen("QEE2.out","a"); fprintf( wfile, "\n[0]\n" ); fclose(wfile);

  do{

      
    U = create_matrix( p+q, 1, EPHEMERAL );
    H = create_matrix( p+q, p+q, EPHEMERAL );

    log_like = 0.0;

wfile=fopen("QEE2.out","a"); fprintf( wfile, "\n start \n" ); fclose(wfile);

    for( c=0; c<nclust; c++ ){

      Ni = Y[c]->nrows;

      if( vlink==1 ){
	gamma = matmult(Z[c],alpha);
	dgda = create_matrix( Ni, 1, EPHEMERAL );
	for( i=0; i<Ni; i++ ){
	  MEL(gamma,i,0) = pow( fabs( MEL(gamma,i,0) ), 0.5 );
	  MEL(dgda,i,0) = 0.5 / MEL(gamma,i,0);
	}
      }
      if( vlink==2 ){
	gamma = matexp( matmult( Z[c], alpha ) );
	dgda = matcopy( gamma );
      }
      make_permanent( gamma );
      make_permanent( dgda );


      Ui = create_matrix(  p+q, 1, EPHEMERAL );
      Hi = create_matrix(  p+q, p+q, EPHEMERAL );

      logPi = 0.0;

      /* 98/06/11 */
      mu_0 = matantilogit( matmult( X[c], beta_0 ) ); 

      /* =-=-=-=-=-=-=-=-=-=-=-=-=-= */
      /*                             */
      /*   Conditional Calcs         */
      /*                             */
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-= */
      if( model==CONDITIONAL ){

	destroy_matrix( Eta[c] );
	Eta[c] = matmult( X[c], beta );
	make_permanent(Eta[c]);

	/* clear marginal mean for recalc */
	for(i=0;i<Ni;i++) MEL( Mu[c], i, 0 ) = 0.0;
	
	for( s=0; s<r; s++ ){

	  logPi_s = 0.0;

	  Di = create_matrix(  p+q, 1, EPHEMERAL );
	  
	  /* evaluate the likelihood kernal for fixed s */
	  for( i=0; i<Ni; i++ ){
	    
	    y_i = MEL( Y[c], i, 0 );
	    n_i = MEL( N[c], i, 0 );
	    
	    eta_i = MEL( Eta[c], i, 0 );
	    
	    gi = MEL( gamma, i, 0 );
	    
	    /* 98/06/11 */
	    /* pm = MEL( mu_0, i, 0 ); */
	    
	    ps = antilogit( eta_i + z[s]*gi );

	    MEL( Mu[c], i, 0 ) += ps * w[s];

	    /* 98/06/11 */	    
	    /* p_q_s = pow( (ps/pm), (double)( y_i ) ) * 
	      pow( (1.-ps)/(1.-pm), (double)(n_i-y_i) ); */

	    p_q_s = y_i * log( ps ) + (n_i-y_i) * log( 1. - ps );
	    
	    logPi_s += p_q_s;

	  }/* i */

	  Pi_s = exp( logPi_s );
	    
	  /* score/hessian contributions for fixed s */
	  for( i=0; i<Ni; i++ ){
	    
	    y_i = MEL( Y[c], i, 0 );
	    n_i = MEL( N[c], i, 0 );
	    
	    eta_i = MEL( Eta[c], i, 0 );
	    
	    gi = MEL( gamma, i, 0 );
	    dgdai = MEL( dgda, i, 0 );
	    
	    ps = antilogit( eta_i + z[s]*gi );

	    for( j=0; j<(p+q); j++ ){
	      
	      if( j < p ){
		dEta_dTheta_j = MEL( X[c], i, j );
	      }else{
		dEta_dTheta_j = MEL( Z[c], i, j-p ) * dgdai * z[s];	  
	      }
	      
	      MEL( Ui, j, 0 ) += Pi_s * (y_i-n_i*ps) * dEta_dTheta_j * w[s];

	      MEL( Di, j, 0 ) += (y_i-n_i*ps) * dEta_dTheta_j;
	      
	      for( k=0; k<(p+q); k++ ){
		
		if( k < p ){
		  dEta_dTheta_k = MEL( X[c], i, k );
		}else{
		  dEta_dTheta_k = MEL( Z[c], i, k-p ) * dgdai * z[s];	  
		}
		
		d2Eta_dTheta_jk = 0.0;
		if( (j>=p) & (k>=p) ){
		  
		  if( vlink==1 ){         /* linear for variance */
		    d2Eta_dTheta_jk = -0.25 * pow( gi, 3.0 )
		      * MEL( Z[c], i, j-p ) 
		      * MEL( Z[c], i, k-p ) * z[s];
		  }
		  if( vlink==2 ){        /* log-linear for std deviation */
		    d2Eta_dTheta_jk = gi * MEL( Z[c], i, j-p ) 
		      * MEL( Z[c], i, k-p ) * z[s];
		  }	
		  
		}/* d2 calc */
		
		MEL( Hi, j, k ) += Pi_s * ( -1.0 * n_i*ps*(1-ps)
					   * dEta_dTheta_j * dEta_dTheta_k
					   + (y_i-n_i*ps) * d2Eta_dTheta_jk ) 
                                   * w[s];
	      }/* k */
	      
	    }/* j */
	    
	  }/* i */

	  logPi += Pi_s * w[s]; 

	  make_permanent( Di );
	  Hi = matadd( Hi, scalar_times_matrix( Pi_s * w[s],
		  matmult( Di, transp(Di) ) ) );
	  destroy_matrix( Di );
	  
	}/* s */
	
      }/* conditional loop */
      
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-= */
      /*                             */
      /*   Marginal Calcs            */
      /*                             */
      /* =-=-=-=-=-=-=-=-=-=-=-=-=-= */
      if( model==MARGINAL ){

	lp = matmult( X[c], beta );
	for( i=0; i<Ni; i++ ) MEL( Mu[c], i, 0 ) = antilogit( MEL(lp,i,0) );
	destroy_matrix( lp );

	/* update deconvolution calcs */

	E = create_matrix( Ni, 5, PERMANENT );

	for( i=0; i<Ni; i++ ){

	  eta_i = MEL( Eta[c], i, 0 );

	  gi = MEL( gamma, i, 0 );

	  mu_i = MEL( Mu[c], i, 0 );
	
	  detadb = 0.0;
	  detads = 0.0;
	  d11 = 0.0;
	  d12 = 0.0;
	  d22 = 0.0;
	
	  flag=0;

	  deconvolve_GH( &eta_i, &detadb, &detads, &d11, &d12, &d22,
		         mu_i, gi, z, w, r, &flag );

	  if( flag==0 ){
	    MEL( Eta[c], i, 0 ) = eta_i;
	  }else{
	    detadb = 0.0;
	    detads = 0.0;
	    d11 = 0.0;
	    d12 = 0.0;
	    d22 = 0.0;	  
	    wfile=fopen("QEE2.out","a");
	    fprintf( wfile, "ERROR:  deconvolve_GH() \n" );
	    fprintf( wfile, "        cluster = %i, subject = %i \n", c, i );
	    fprintf( wfile, "        mu = %f, eta = %f\n", mu_i, eta_i );
	    fprintf( wfile, "        derivatives set to 0 for this cycle\n");
	    fclose(wfile);
	  }

	  MEL( E, i, 0 ) = detadb;
	  MEL( E, i, 1 ) = detads;
	  MEL( E, i, 2 ) = d11;
	  MEL( E, i, 3 ) = d12;
	  MEL( E, i, 4 ) = d22; 

	}/* i */

	for( s=0; s<r; s++ ){

	  logPi_s = 0.0;

	  Di = create_matrix( p+q, 1, EPHEMERAL );

	  /* evaluate the likelihood kernal for fixed s */
	  for( i=0; i<Ni; i++ ){
	    
	    y_i = MEL( Y[c], i, 0 );
	    n_i = MEL( N[c], i, 0 );
	    
	    eta_i = MEL( Eta[c], i, 0 );
	    
	    gi = MEL( gamma, i, 0 );
	    
	    /* 98/06/11 */
	    /* pm = MEL( mu_0, i, 0 ); */
	    
	    ps = antilogit( eta_i + z[s]*gi );

	    /* 98/06/11 */
	    /* p_q_s = pow( (ps/pm), (double)( y_i ) ) * 
	      pow( (1.-ps)/(1.-pm), (double)(n_i-y_i) ); */

	    p_q_s = y_i * log( ps ) + (n_i-y_i) * log( 1. - ps );
	    
	    logPi_s += p_q_s;

	  }/* i */

	  Pi_s = exp( logPi_s );

	  /* score / hessian calculations */
	  for( i=0; i<Ni; i++ ){

	    y_i = MEL( Y[c], i, 0 );
	    n_i = MEL( N[c], i, 0 );
	    
	    eta_i = MEL( Eta[c], i, 0 );
	    
	    gi = MEL( gamma, i, 0 );
	    dgdai = MEL( dgda, i, 0 );
	    
	    /* 98/06/11 */
	    /* pm = MEL( mu_0, i, 0 ); */

	    ps = antilogit( eta_i + z[s]*gi );

	    detadb = MEL( E, i, 0 );
	    detads = MEL( E, i, 1 );
	    d11 = MEL( E, i, 2 );
	    d12 = MEL( E, i, 3 );
	    d22 = MEL( E, i, 4 );
	  
	    for( j=0; j<(p+q); j++ ){
	    
	      if( j < p ){
		dEta_dTheta_j = detadb * MEL( X[c], i, j );
	      }else{
		dEta_dTheta_j = (detads + z[s]) * MEL( Z[c], i, j-p ) * dgdai;
	      }
	    
	      MEL( Ui, j, 0 ) += Pi_s * (y_i-n_i*ps) * dEta_dTheta_j * w[s];

	      MEL( Di, j, 0 ) += (y_i-n_i*ps)* dEta_dTheta_j;
	     
	      for( k=0; k<(p+q); k++ ){
	      
		if( k < p ){
		  dEta_dTheta_k = detadb * MEL( X[c], i, k );
		}else{
		  dEta_dTheta_k = (detads + z[s]) * MEL( Z[c], i, k-p ) * 
		                     dgdai;
		}
	      
		d2Eta_dTheta_jk = 0.0;

		if( (j<p) & (k<p) ){
		
		  d2Eta_dTheta_jk = d11 * MEL( X[c], i, j ) * 
		                          MEL( X[c], i, k );
		
		}/* d2 calc */

		if( (j<p) & (k>=p) ){
		
		  d2Eta_dTheta_jk = d12 * MEL( X[c], i, j ) * 
		                          MEL( Z[c], i, k-p) * dgdai;
		}/* d2 calc */
		
		if( (j>=p) & (k<p) ){
		  
		  d2Eta_dTheta_jk = d12 * MEL( X[c], i, k ) * 
		                          MEL( Z[c], i, j-p) * dgdai;
		}/* d2 calc */

		if( (j>=p) & (k>=p) ){
		
		  if( vlink==1 ){         /* linear for variance */
		    d2Eta_dTheta_jk = -0.25 * pow( gi, 3.0 )
		      * MEL( Z[c], i, j-p ) 
		      * MEL( Z[c], i, k-p ) * ( z[s]  + detads )
		      + d22 * MEL( Z[c], i, j-p ) * MEL( Z[c], i, k-p ) *
		        pow( dgdai, 2.0 );
		  }
		  if( vlink==2 ){        /* log-linear for std deviation */
		    d2Eta_dTheta_jk = gi * MEL( Z[c], i, j-p ) 
		      * MEL( Z[c], i, k-p ) * ( z[s]  + detads )
		      + d22 * MEL( Z[c], i, j-p )*MEL( Z[c], i, k-p ) *
		        pow( dgdai, 2.0 );
		  }	
		
		}/* d2 calc */
	      
		MEL( Hi, j, k ) += Pi_s * ( -1.0 * n_i*ps*(1-ps) 
					   * dEta_dTheta_j * dEta_dTheta_k
					   + (y_i-n_i*ps) * d2Eta_dTheta_jk ) 
		                           * w[s];
	      }/* k */
	    
	    }/* j */

	  }/* i */

	  logPi += Pi_s * w[s];

	  make_permanent( Di );
	  Hi = matadd( Hi, scalar_times_matrix( Pi_s * w[s],
		  matmult( Di, transp(Di) ) ) );
	  destroy_matrix( Di );

      	}/* s */

	destroy_matrix( E ); 

      }/* marginal loop */

      for( j=0; j<(p+q); j++ ){
	MEL( Ui, j, 0 ) = MEL(Ui,j,0)/logPi;
      }

      for( j=0; j<(p+q); j++ ){
	for( k=0; k<(p+q); k++ ){

	  if( iter >= NSAFE ){
	    MEL( Hi, j, k ) = MEL( Ui, j, 0 ) * MEL( Ui, k, 0 ) -
	                    MEL( Hi, j, k )/logPi; 
	  }else{
	    MEL( Hi, j, k ) = MEL( Ui, j, 0 ) * MEL( Ui, k, 0 ); 
	  }

	}
      }

      U = matadd( Ui, U );
      H = matadd( Hi, H );

      /* 98/06/11 */
      /*
      logPi_indep = 0.0;
      for(i=0;i<Ni;i++){
	pm = MEL( mu_0, i, 0 );
	logPi_indep += MEL(Y[c],i,0)*log( pm/(1.-pm) );
	logPi_indep += MEL(N[c],i,0)*log(1.-pm);
      }

      log_like += log( logPi ) + logPi_indep;
      */
      log_like += log( logPi );

      destroy_matrix( mu_0 );
      destroy_matrix( gamma );
      destroy_matrix( dgda );

    }/* c */
	      

wfile=fopen("QEE2.out","a"); fprintf( wfile, "\n (2) \n" ); fclose(wfile);

/* =-=-=-=-=-=-=-=-= */
/*   regularization  */
/* =-=-=-=-=-=-=-=-= */

    for( k=0; k<p+q; k++ ){

      if( k<p ){
	MEL( U, k, 0 ) -= lambda * MEL( beta, k, 0 );	
      }else{
	MEL( U, k, 0 ) -= lambda * MEL( alpha, k-p, 0 );	
      }

      MEL( H, k, k ) += lambda + RIDGE;

    }

    step=1.0;
    /*    if( HALFSTEP ) step = 1.0 - pow( 0.75, (double)(iter+1) ); */
    if( HALFSTEP ) step = 0.5;

    make_permanent( U );
    make_permanent( H );

    wfile=fopen("QEE2.out","a");
    fprintf(wfile,"\n U \n");
    fmatdump(wfile,transp(U));
    fprintf(wfile,"\n H \n");
    fmatdump(wfile,transp(H));
    fclose(wfile);

    delta = matmult( luinv(H), U );
    for( i=0; i<p; i++ ) MEL(beta,i,0) += step*MEL(delta,i,0);
    for( i=0; i<q; i++ ) MEL(alpha,i,0) += step*MEL(delta,p+i,0);
    dmax = matmaxabs( delta );

    iter++;

    if( dmax < tolerance ) converge=1;
	  
    wfile=fopen("QEE2.out","a");
    fprintf(wfile,"\n beta \n");
    fmatdump(wfile,transp(beta));
    fprintf(wfile,"\n alpha \n");
    fmatdump(wfile,transp(alpha));
    fprintf(wfile,"\n dmax=%f, logL=%f\n",dmax,log_like);
    fclose(wfile);

    /* =-=-=-=-=-=-=- */
    /*  cleaning...   */
    /* =-=-=-=-=-=-=- */

    if( !converge ){
    /*  destroy_matrix( gamma );
      destroy_matrix( dgda );*/
      destroy_matrix( U );
      destroy_matrix( H );
    }
wfile=fopen("QEE2.out","a"); fprintf( wfile, "\n end \n" ); fclose(wfile);

  }while( !converge && (iter < maxiter) );

wfile=fopen("QEE2.out","a"); fprintf( wfile, "\n[1]\n" ); fclose(wfile);


  /* add the appropriate cleaning... */


  to_S( beta, S_beta );

  to_S( alpha, S_alpha );

  if( model==CONDITIONAL ){
    count = 0;
    for(c=0;c<nclust;c++){
      Ni = Mu[c]->nrows;
      for(i=0;i<Ni;i++){
	*(S_eta + count ) = logit( MEL( Mu[c], i, 0 ) );
	count++;
      }
    }
  }
  if( model==MARGINAL ){
    count = 0;
    for(c=0;c<nclust;c++){
      Ni = Mu[c]->nrows;
      for(i=0;i<Ni;i++){
	*(S_eta + count ) = MEL( Eta[c], i, 0 );
	count++;
      }
    }
  }

  if( converge ){
    for( k=0; k<(p+q); k++ ) MEL( H, k, k ) -= RIDGE;
    H = luinv( H );
    to_S( H, S_modcov );
    *S_logL = log_like;
  }

  *S_tol = dmax;

  *S_maxiter = iter;

}/* end of routine */

/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
int split( matptr, discptr , matarrptr )
MATRIX *matptr, *discptr, *matarrptr[];
{   /* discriminator vector assumed to be integer-valued dbls */
int i, j, istart, k, start, end;
if ( discptr->ncols != 1 )
	{
	fprintf(stderr,"split: discriminator must be column vec.\n");
	fprintf(stderr,"split: ncols = %d.\n", discptr->ncols);
	fprintf(stderr,"split: fatal error.\n");
	Seterr_and_terminate(SPLIT_FAIL );
	}

k = 0;

istart = (int)MEL( discptr , 0 , 0 );
start = 0;
end = 0;
for ( i = 1 ; i <= discptr->nrows ; i++ )
	{
	if (( MEL( discptr , i, 0 ) != istart ) ||
			i == (discptr->nrows ) )
		{
		matarrptr[k] = matcopy( extract_rows( matptr , start, end ) );
		make_permanent( matarrptr[k] );
		k++;
		start = end+1;
		istart = MEL( discptr, i, 0 );
		}
	if (start < discptr->nrows ) end++ ;
	}
/* DOES NOT CLEAN */
return k;
}
 
MATRIX *stack( x, y )
MATRIX *x, *y;
{
int xc, yc, xr, yr, i, j;
MATRIX *tmp;

xc = x->ncols;
yc = y->ncols;
if ( xc != yc ) 
	{
	fprintf(stderr, "M+-: stack: incompatible columns.\n");
	Seterr_and_terminate( CAN_T_STACK_MATRICES )
	}
xr = x->nrows;
yr = y->nrows;
tmp = create_matrix( xr+yr, xc, EPHEMERAL );
for ( i = 0 ; i < xr+yr ; i++ )
	{
	for ( j = 0 ; j < xc ; j++ )
		{
		MEL(tmp,i,j) = ( i >= xr ) ? MEL( y, i-xr , j ) : MEL( x, i, j );
		}
	}
free_if_ephemeral(x);
free_if_ephemeral(y);
return tmp;
}


double *get_GH_z( q )

  int q;
{
  double *z;

  if( (q!=3)*(q!=5)*(q!=10)*(q!=20)*(q!=50) ){
    printf("get_GH_z():  Error.  q should be 3, 5, 10, 20, or 50.\n");
    exit( 0 );
  }

  switch( q ){
  case 3:
    z = (double *)calloc( 3, (unsigned)sizeof( double ) );
    z[0] = -1.224744871391589 * sqrt(2.0);
    z[1] =  0.0;
    z[2] =  1.224744871391589 * sqrt(2.0);
  case 5:
    z = (double *)calloc( 5, (unsigned)sizeof( double ) );
    z[0] = -2.020182870456086 * sqrt(2.0);
    z[1] = -0.958572464613819 * sqrt(2.0);
    z[2] = 0.000000000000000 * sqrt(2.0);
    z[3] = 0.958572464613819 * sqrt(2.0);
    z[4] = 2.020182870456086 * sqrt(2.0);
    break;
  case 10:
    z = (double *)calloc( 10, (unsigned)sizeof( double ) );
    z[0] = -3.436159118837738 * sqrt(2.0);
    z[1] = -2.532731674232790 * sqrt(2.0);
    z[2] = -1.756683649299882 * sqrt(2.0);
    z[3] = -1.036610829789514 * sqrt(2.0);
    z[4] = -0.342901327223705 * sqrt(2.0);
    z[5] = 0.342901327223705 * sqrt(2.0);
    z[6] = 1.036610829789514 * sqrt(2.0);
    z[7] = 1.756683649299882 * sqrt(2.0);
    z[8] = 2.532731674232790 * sqrt(2.0);
    z[9] = 3.436159118837738 * sqrt(2.0);
    break;
  case 20:
    z = (double *)calloc( 20, (unsigned)sizeof( double ) );
    z[0] = -5.3874808900112 * sqrt(2.0);
    z[1] = -4.6036824495507 * sqrt(2.0);
    z[2] = -3.9447640401156 * sqrt(2.0);
    z[3] = -3.3478545673832 * sqrt(2.0);
    z[4] = -2.7888060584281 * sqrt(2.0);
    z[5] = -2.2549740020892 * sqrt(2.0);
    z[6] = -1.7385377121166 * sqrt(2.0);
    z[7] = -1.2340762153953 * sqrt(2.0);
    z[8] = -0.7374737285454 * sqrt(2.0);
    z[9] = -0.2453407083009 * sqrt(2.0);
    z[10] = 0.2453407083009 * sqrt(2.0);
    z[11] = 0.7374737285454 * sqrt(2.0);
    z[12] = 1.2340762153953 * sqrt(2.0);
    z[13] = 1.7385377121166 * sqrt(2.0);
    z[14] = 2.2549740020892 * sqrt(2.0);
    z[15] = 2.7888060584281 * sqrt(2.0);
    z[16] = 3.3478545673832 * sqrt(2.0);
    z[17] = 3.9447640401156 * sqrt(2.0);
    z[18] = 4.6036824495507 * sqrt(2.0);
    z[19] = 5.3874808900112 * sqrt(2.0);
    break;
  case 50:
    z = (double *)calloc( 50, (unsigned)sizeof( double ) );
    z[0] = -12.9858845;
    z[1] = -12.0530184;
    z[2] = -11.2792333;
    z[3] = -10.5873817;
    z[4] = -9.9480357;
    z[5] = -9.3460396;
    z[6] = -8.7722996;
    z[7] = -8.2208159;
    z[8] = -7.6873624;
    z[9] = -7.1688148;
    z[10] = -6.6627754;
    z[11] = -6.1673474;
    z[12] = -5.6809923;
    z[13] = -5.2024350;
    z[14] = -4.7305986;
    z[15] = -4.2645578;
    z[16] = -3.8035057;
    z[17] = -3.3467278;
    z[18] = -2.8935827;
    z[19] = -2.4434875;
    z[20] = -1.9959047;
    z[21] = -1.5503332;
    z[22] = -1.1062993;
    z[23] = -0.6633497;
    z[24] = -0.2210452;
    z[25] = 0.2210452;
    z[26] = 0.6633497;
    z[27] = 1.1062993;
    z[28] = 1.5503332;
    z[29] = 1.9959047;
    z[30] = 2.4434875;
    z[31] = 2.8935827;
    z[32] = 3.3467278;
    z[33] = 3.8035057;
    z[34] = 4.2645578;
    z[35] = 4.7305986;
    z[36] = 5.2024350;
    z[37] = 5.6809923;
    z[38] = 6.1673474;
    z[39] = 6.6627754;
    z[40] = 7.1688148;
    z[41] = 7.6873624;
    z[42] = 8.2208159;
    z[43] = 8.7722996;
    z[44] = 9.3460396;
    z[45] = 9.9480357;
    z[46] = 10.5873817;
    z[47] = 11.2792333;
    z[48] = 12.0530184;
    z[49] = 12.9858845;
    break;
  default:
    printf("get_GH_z():  Error q unknown.\n");
    exit( 0 );
  }

  return z;
}
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

double *get_GH_w( q )

  int  q;
{
  double *w;

  if( (q!=3)*(q!=5)*(q!=10)*(q!=20)*(q!=50) ){
    printf("get_GH_z():  Error.  q should be 3, 5, 10, 20, or 50.\n");
    exit( 0 );
  }

  switch( q ){
  case 3:
    w = (double *)calloc( 3, (unsigned)sizeof( double ) );
    w[0] = ( 2.954089751509e-1 ) /sqrt(PI);
    w[1] = ( 1.181635900604 ) /sqrt(PI);
    w[2] = ( 2.954089751509e-1 ) /sqrt(PI);
  case 5:
    w = (double *)calloc( 5, (unsigned)sizeof( double ) );
    w[0] = 1.995324205905e-2 / sqrt(PI);
    w[1] = 3.936193231522e-1 / sqrt(PI);
    w[2] = 9.453087204829e-1 / sqrt(PI);
    w[3] = 3.936193231522e-1 / sqrt(PI);
    w[4] = 1.995324205905e-2 / sqrt(PI);
    break;
  case 10:
    w = (double *)calloc( 10, (unsigned)sizeof( double ) );
    w[0] =   7.640432855233e-6 / sqrt(PI);
    w[1] =   1.343645746781e-3 / sqrt(PI);
    w[2] =   3.387439445548e-2 / sqrt(PI);
    w[3] =   2.401386110823e-1 / sqrt(PI);
    w[4] =   6.108626337353e-1 / sqrt(PI);
    w[5] =  6.108626337353e-1 / sqrt(PI);
    w[6] =  2.401386110823e-1 / sqrt(PI);
    w[7] =  3.387439445548e-2 / sqrt(PI);
    w[8] =  1.343645746781e-3 / sqrt(PI);
    w[9] =  7.640432855233e-6 / sqrt(PI);
    break;
  case 20:
    w = (double *)calloc( 20, (unsigned)sizeof( double ) );
    w[0] =  2.229393645534e-13 /sqrt(PI);
    w[1] =  4.399340992273e-10 /sqrt(PI);
    w[2] =  4.399340992273e-10 /sqrt(PI);
    w[3] =  7.802556478532e-6  /sqrt(PI);
    w[4] =  2.283386360163e-4  /sqrt(PI);
    w[5] =  3.243773342238e-3  /sqrt(PI);
    w[6] =  2.481052088746e-2  /sqrt(PI);
    w[7] =  1.090172060200e-1  /sqrt(PI);
    w[8] =  2.866755053628e-1  /sqrt(PI);
    w[9] =  4.622436696006e-1  /sqrt(PI);
    w[10] = 4.622436696006e-1  /sqrt(PI);
    w[11] = 2.866755053628e-1  /sqrt(PI);
    w[12] = 1.090172060200e-1  /sqrt(PI);
    w[13] = 2.481052088746e-2  /sqrt(PI);
    w[14] = 3.243773342238e-3  /sqrt(PI);
    w[15] = 2.283386360163e-4  /sqrt(PI);
    w[16] = 7.802556478532e-6  /sqrt(PI);
    w[17] = 1.086069370769e-7  /sqrt(PI);
    w[18] = 4.399340992273e-10 /sqrt(PI);
    w[19] = 2.229393645534e-13 /sqrt(PI);
    break;
  case 50:
    w = (double *)calloc( 50, (unsigned)sizeof( double ) );
    w[0] = 1.034608e-37;
    w[1] = 9.443415e-33;
    w[2] = 6.856281e-29;
    w[3] = 1.206045e-25;
    w[4] = 7.995094e-23;
    w[5] = 2.522483e-20;
    w[6] = 4.368172e-18;
    w[7] = 4.566698e-16;
    w[8] = 3.083829e-14;
    w[9] = 1.414229e-12;
    w[10] = 4.576637e-11;
    w[11] = 1.077061e-09;
    w[12] = 1.888226e-08;
    w[13] = 2.514610e-07;
    w[14] = 2.584938e-06;
    w[15] = 2.078485e-05;
    w[16] = 1.321726e-04;
    w[17] = 6.708281e-04;
    w[18] = 2.738161e-03;
    w[19] = 9.045054e-03;
    w[20] = 2.430481e-02;
    w[21] = 5.334352e-02;
    w[22] = 9.593054e-02;
    w[23] = 1.416854e-01;
    w[24] = 1.721259e-01;
    w[25] = 1.721259e-01;
    w[26] = 1.416854e-01;
    w[27] = 9.593054e-02;
    w[28] = 5.334352e-02;
    w[29] = 2.430481e-02;
    w[30] = 9.045054e-03;
    w[31] = 2.738161e-03;
    w[32] = 6.708281e-04;
    w[33] = 1.321726e-04;
    w[34] = 2.078485e-05;
    w[35] = 2.584938e-06;
    w[36] = 2.514610e-07;
    w[37] = 1.888226e-08;
    w[38] = 1.077061e-09;
    w[39] = 4.576637e-11;
    w[40] = 1.414229e-12;
    w[41] = 3.083829e-14;
    w[42] = 4.566698e-16;
    w[43] = 4.368172e-18;
    w[44] = 2.522483e-20;
    w[45] = 7.995094e-23;
    w[46] = 1.206045e-25;
    w[47] = 6.856281e-29;
    w[48] = 9.443415e-33;
    w[49] = 1.034608e-37;
    break;
  default:
    printf("get_GH_w():  Error.  q unknown.\n");
    exit( 0 );
  }

  return w;
}
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

double antilogit( x )

  double x;
{
  double out;

  out = exp( x )/(1.0+exp( x ) );

  return out;

}
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

double logit( x )

  double x;
{
  double out;

  out = log( x ) - log( 1.-x );

  return out;

}
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

void deconvolve_GH( eta, detadb, detads, d2etadb2, d2etadbds, d2etads2,
		    mu, sigma, z, w, r, flag)

  double  *eta, *detadb, *detads, *d2etadb2, *d2etadbds, *d2etads2;
  double  mu, sigma, *z, *w;
  integer *flag, r;

{

  double   etai, dmuidetai, dmuids, new_mui, delmu, ps;
  double   d1, d2, int2, int3, int4;
  integer  eta_converge, s, count;       
  etai = *eta;

  eta_converge = 0;
  count = 0;

  do{
    dmuidetai = 0.0;
    dmuids = 0.0;
    new_mui = 0.0;
    int2 = 0.0;
    int3 = 0.0;
    int4 = 0.0;
    for( s=0; s<r; s++ ){
      ps = antilogit( etai + sigma*z[s] );
      new_mui += w[s] * ps;
      dmuidetai += w[s] * ps * (1.-ps);
      dmuids += w[s] * ps * (1.-ps) * z[s];
      int2 += w[s] * (1.-2.*ps)*ps*(1.-ps);
      int3 += w[s] * (1.-2.*ps)*ps*(1.-ps) * z[s];
      int4 += w[s] * (1.-2.*ps)*ps*(1.-ps) * z[s] * z[s];
    }
    delmu = ( new_mui - mu )/( dmuidetai );
    etai -= delmu;
    if( fabs(delmu) < ETA_TOLERANCE ) eta_converge = 1;
    count++;
  }while( count < ETA_MAX_ITER && !eta_converge );
  
  if( count >= ETA_MAX_ITER ){
    *flag = 1;
  }else{
    *flag = 0;
  }

  *eta = etai;

  if( dmuidetai < ETA_EPS ) dmuidetai = ETA_EPS;

  d1 = ( mu*(1.0-mu) )/dmuidetai;

  d2 = -1.0 * dmuids/dmuidetai;

  *detadb = d1;

  *detads = d2;

  *d2etadb2 = ( mu*(1.-mu)*(1.-2.*mu) - d1*d1*int2 )/dmuidetai;

  *d2etadbds = -1.0 * ( d1*d2*int2 + d1*int3 )/dmuidetai;

  *d2etads2 = -1.0 * ( d2*d2*int2 + 2.*d2*int3 + int4 )/dmuidetai; 

  return;

}

