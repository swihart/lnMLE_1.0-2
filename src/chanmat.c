
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
/* This point denotes the beginning of file chanmat.c :: B.A.C. */
#include "chanmatstruct.h"
int choose2(n) int n;
{
return (n*(n-1)/2);
}
MATRIX *create_matrix( nrows, ncols, permanence )
/* $Y = |create_matrix|(r,c,\cdot) :\Rightarrow Y \in M_{r\times c} \wedge
Y = 0 $ */
int nrows, ncols, permanence;
{
MATRIX *tmp;
double *head;
int i;

tmp = (MATRIX *) calloc ( 1, sizeof ( struct matrix ) );

if ( tmp == NULL )
	{
	fprintf( stderr , "create_matrix: malloc attempt %d d.\n",
				sizeof( struct matrix ));
	Seterr_and_terminate( NO_MEM_MATSTRUCT );
	}

tmp->nrows = nrows;
tmp->ncols = ncols;
tmp->permanence = permanence;

tmp->data = ( double * ) calloc ( 1,  nrows * ncols * sizeof ( double ) ) ;

if ( tmp->data == NULL )
	{
	fprintf( stderr , "create_matrix: malloc attempt %d d.\n",
				(unsigned)nrows*ncols);
	fprintf( stderr , "create_matrix: nrows=%d ncols=%d.\n",
				nrows, ncols );
	Seterr_and_terminate( NO_MEM_MATDATA );
	}

head = tmp->data;
for ( i = 0 ; i < nrows*ncols ; i++ )
	{
	*head = (double)0.;
	head++;
	}

return tmp;
}

void destroy_matrix( mat )
MATRIX *mat;
{
if (mat == (MATRIX *) NULL ) return;
mat->nrows = 0;
mat->ncols = 0;
if (mat->data != (double *)NULL ) free((char *) mat->data );
mat->data = (double *)NULL;
if (mat != (MATRIX *)NULL) free((char *) mat );
mat = (MATRIX *)NULL;
}

 

MATRIX *transp( mat )
/*
$Y = transp(X_{r\times c}) :\Rightarrow Y \in M_{c\times r} \wedge
y_{ji} = x_{ij}. $
*/
MATRIX *mat;
{
double *telem, *inelem, *tbase;
int nelem;
MATRIX *tmp;

tmp = create_matrix( mat->ncols, mat->nrows , EPHEMERAL );
inelem = mat->data;
tbase = tmp->data;
telem = tbase;
for ( nelem = 0 ; nelem < ( mat->ncols * mat->nrows ) ; nelem++ )
	{
	*telem = *(inelem++);
	telem += mat->nrows;
	if ( nelem % mat->ncols == (mat->ncols)-1 )
		telem = ++tbase;
	}
if ( is_ephemeral( mat ) ) destroy_matrix( mat );
return tmp;

}
 
MATRIX *corner( mat, nr, nc )
/*
$ r \leq r^{\prime} \wedge c \leq c^{\prime} 
\wedge X \in M_{r^{\prime} \times c^{\prime}} 
\wedge Y = $|corner(X,r,c)|$ :\Rightarrow $ */
/* $ Y \in M_{r \times c} \wedge y_{ij} = x_{ij}, 
\thinspace 1 \leq i \leq r, \thinspace 1 \leq j \leq c $ */

MATRIX *mat;
int nr, nc;
{
MATRIX *tmp;
double *load;
int i,j,sr, sc;
sr = mat->nrows;
sc = mat->ncols;
if ((nr > sr) || (nc > sc))
	{
	fprintf( stderr, "corner: request not a submatrix.\n");
	fprintf( stderr, "corner: fatal error.\n");
	Seterr_and_terminate(CORNER_FAIL);
	}
tmp = create_matrix( nr, nc, EPHEMERAL );
load = tmp->data;
for ( i = 0 ; i < nr ; i++ )
	{
	for ( j = 0 ; j < nc ; j++ )
		{
		*(load++) = MEL(mat, i, j);
		}
	}
free_if_ephemeral( mat );
return tmp;
}


MATRIX *extract_rows(Source,start,end)
	/* purely zero-based */
	
MATRIX *Source;
int start, end;
{
MATRIX *temp;
int rows_to_get, i, j;

if( Source->nrows < end ){
  printf("chanmat: extract_rows() end is beyond rows of Source!\n");
  exit( 0 ); }

rows_to_get = end - start + 1;

if( rows_to_get <= 0 ){
  printf("chanmat: extract_rows() rows_to_get = %i\n",rows_to_get);
  exit( 0 ); }

temp = create_matrix(rows_to_get,Source->ncols,EPHEMERAL);

for ( i = 0 ; i < rows_to_get ; i++ )
	{
	for ( j = 0 ; j < Source->ncols ; j++ )
		{
		*(ELREF(temp,i,j)) = *(ELREF(Source,start,j));
		}
	start++;
	}
/* DOES NOT CLEAN */
return temp;
}

MATRIX *extract_cols( x , start , end )
MATRIX *x;
int start, end;
{
MATRIX *tmp;
tmp = transp(x);
tmp = extract_rows( tmp, start, end );
tmp = transp(tmp);
free_if_ephemeral(x);
return tmp;
}

MATRIX *matcopy(inmat)
MATRIX *inmat; 
{
int i, j;
MATRIX *outmat;

outmat = create_matrix(inmat->nrows,inmat->ncols,EPHEMERAL);
for ( i = 0 ; i < inmat->nrows ; i++ )
	{
	for ( j = 0 ; j < inmat->ncols ; j++ )
		{
		*(ELREF(outmat,i,j)) = *(ELREF(inmat,i,j));
		}
	}
/* DOES NOT CLEAN */
return outmat;
}
 
void plug( plugm, socket, row, col )
int row, col;
MATRIX *plugm, *socket;  /* not a unix socket */
{
int pcol, prow;
double *sockload, *plughead, *sockrow_start;
int i,j;

pcol = plugm->ncols;
prow = plugm->nrows;

if ( pcol+col > socket->ncols || prow+row > socket->nrows )
	{
	fprintf( stderr,"M+-: plug: socket too small. Dies.\n");
	Seterr_and_terminate(PLUG_FAIL);
	}

sockload = socket->data + col + row*(socket->ncols);
plughead = plugm->data;
sockrow_start = sockload;

for ( i = 0 ; i < prow ; i++ )
	{
	sockload = sockrow_start;
	for ( j = 0 ; j < pcol ; j++ )
		{
		*(sockload++) = *(plughead++ );
		}
	sockrow_start += socket->ncols;
	}
free_if_ephemeral( plugm );
}
 
MATRIX *form_diag( vec )
MATRIX *vec;
{
MATRIX *tmp;
int i, ord;

ord = vec->nrows;
tmp = create_matrix( ord,  ord, EPHEMERAL );
for ( i = 0 ; i < ord ; i++ )
	*(ELREF(tmp,i,i)) = MEL(vec,i,0);
free_if_ephemeral( vec );
return tmp;
}

MATRIX *band( in, wid )
MATRIX *in;
int wid;
{
MATRIX *tmp;
int i, j;
tmp = matcopy( in );
for ( i = 0 ; i < in->nrows ; i++ )
	{
	for ( j = i+wid ; j < in->ncols ; j++ )
		{
		MEL( tmp, i, j ) = (double)0.;
		if (( i < in->ncols ) && ( j < in->nrows ))
			{
			MEL( tmp, j, i ) = (double)0.;
			}
		}
	}
free_if_ephemeral( in );
return tmp;
}

MATRIX *toeplitz( in )
MATRIX *in;
{
MATRIX *toep, *tin, *tmp;
int n, p, inrows, incols, i, j;

inrows = in->nrows;
incols = in->ncols;

if ( ( inrows > incols ) ? inrows % incols : incols % inrows )
	{
	fprintf(stderr,"M+-:toeplitz: argument invalid. Dies.\n");
	Seterr_and_terminate(BAD_TOEPLITZ_ARG);
	}

if ( inrows > incols )
	{
	p = incols;
	n = inrows/p;
	tin = matcopy(in);
	free_if_ephemeral(in);
	}
else
	{
	p = inrows;
	n = incols/p;
	tin = transp(in);
	}

toep = create_matrix( n*p, n*p, EPHEMERAL );

for ( i = 0 ; i < n ; i ++ )
	{
	tmp = extract_rows( tin, i*p, (i*p)+p-1 );
	make_permanent(tmp);
	if ( i == 0 )
		{
		for ( j = 0 ; j < n ; j++ )
			{
			if ( inrows > incols )
				plug( tmp, toep, j*p, j*p );
			else
				plug( transp(tmp), toep, j*p, j*p );
			}
		}
	else
		{
		for ( j = 0 ; j < n-i ; j++ )
			{
			plug( transp(tmp), toep, j*p, (j+i)*p );
			plug( tmp, toep, (j+i)*p, j*p );
			}
		}
	destroy_matrix(tmp);
	}
destroy_matrix(tin);
return toep;
}

MATRIX *star( x ) MATRIX *x;
{
MATRIX *tmp;
int xr, xc, nc2, i, j, curr;
xc = x->ncols;
if ( xc > 1 ) fprintf(stderr,"M+-: star: should have colvec.\n");
xr = x->nrows;
nc2 = choose2(xr);

tmp = create_matrix( nc2, 1, EPHEMERAL );

curr = 0;
for ( i = 0 ; i < (xr-1) ; i++ )
	{
	for ( j = (i+1); j < xr ; j++ )
		{
		MEL( tmp, curr, 0 ) = MEL( x, i, 0 );
		curr++;
		}
	}
free_if_ephemeral( x );
return tmp;
}
MATRIX *tilde( x ) MATRIX *x;
{
MATRIX *tmp;
int xr, xc, nc2, i, j, curr;
xc = x->ncols;
if ( xc > 1 ) fprintf(stderr,"M+-: tilde: should have colvec.\n");
xr = x->nrows;
nc2 = choose2(xr);

tmp = create_matrix( nc2, 1, EPHEMERAL );

curr = 0;
for ( i = 0 ; i < (xr-1) ; i++ )
	{
	for ( j = (i+1); j < xr ; j++ )
		{
		MEL( tmp, curr, 0 ) = MEL( x, j, 0 );
		curr++;
		}
	}
free_if_ephemeral( x );
return tmp;
}

 
 

#define get_nelem( x ) (((x)->nrows) * ((x)->ncols))

double elsum( x )
MATRIX *x;
{
double t=0.;
double *loc;
int i, nelem;

nelem = get_nelem( x );
loc = x->data;
for ( i = 0 ; i < nelem ; i++ )
	t += *(loc++);
if ( is_ephemeral(x) ) destroy_matrix(x);
return t;
}

MATRIX *matabs( x )
MATRIX *x;
{
double *load, *look;
double fabs();
MATRIX *tmp;
int nelem, i;

nelem = get_nelem( x );
tmp = create_matrix( x->nrows, x->ncols , EPHEMERAL );
load = tmp->data;
look = x->data;
for ( i = 0 ; i < nelem ; i++ )
	*(load++) = fabs(*look++);
free_if_ephemeral(x);
return tmp ;
}

double matmax( x )
MATRIX *x;
{
double t;
double *loc;
int i, nelem;

nelem = get_nelem( x );
loc = x->data;
t = MEL(x,0,0);
for ( i = 0 ; i < nelem ; i++ )
	{
	if ( *(loc) > t ) t = *(loc);
	loc++;
	}
free_if_ephemeral( x );
return t;
}

double matmaxabs( x ) MATRIX* x;
{
double tmp;
int xr, xc, i, j;
xr = x->nrows;
xc = x->ncols;
tmp = fabs(MEL(x,0,0));
for ( i = 0 ; i < xr; i++ )
	{
	for ( j = 0 ; j < xc ; j++ )
		{
		if ( fabs(MEL(x,i,j)) > tmp ) tmp = fabs(MEL(x,i,j));
		}
	}
/* cleaning added: pjh 95/07/16 */
free_if_ephemeral( x );

return tmp;
}

MATRIX *matexp( x )
MATRIX *x;
{
double *load, *look;
double exp();
MATRIX *tmp;
int nelem, i;

nelem = get_nelem( x );
tmp = create_matrix( x->nrows, x->ncols , EPHEMERAL );
load = tmp->data;
look = x->data;
for ( i = 0 ; i < nelem ; i++ )
	*(load++) = exp(*look++);
free_if_ephemeral(x);
return tmp ;
}

MATRIX *matantilogit( x )
MATRIX *x;
{
double *load, *look;
double exp();
MATRIX *tmp;
int nelem, i;

nelem = get_nelem( x );
tmp = create_matrix( x->nrows, x->ncols , EPHEMERAL );
load = tmp->data;
look = x->data;
for ( i = 0 ; i < nelem ; i++ )
	{
	*(load++) = exp(*look)/(1.+exp(*look));
	look++;
	}
free_if_ephemeral(x);
return tmp ;
}

MATRIX *oneminus( x )
MATRIX *x;
{
double *load, *look;
double exp();
MATRIX *tmp;
int nelem, i;

nelem = get_nelem( x );
tmp = create_matrix( x->nrows, x->ncols , EPHEMERAL );
load = tmp->data;
look = x->data;
for ( i = 0 ; i < nelem ; i++ )
	*(load++) = 1.-(*look++);
free_if_ephemeral(x);
return tmp ;
}
 

MATRIX *matadd( mat1, mat2 )
MATRIX *mat1, *mat2;
{
MATRIX *result;
double *mat1base, *mat2base, *resbase;
int i, j, nlen,z=0;
if ( ( mat1->ncols != mat2->ncols ) || ( mat1->nrows != mat2->nrows ) )
	{
	fprintf(stderr,"matadd: args (%dx%d) + (%dx%d) don't conform.\n",
		mat1->nrows, mat1->ncols, mat2->nrows,
		mat2->ncols );
	fprintf(stderr,"matadd: fatal error.  exits. \n");
	Seterr_and_terminate(MATADD_NONCONFORMITY);
	}
result = create_matrix( mat1->nrows , mat1->ncols, EPHEMERAL);
resbase = result->data;
mat1base = mat1->data;
mat2base = mat2->data;
for ( j = 0 ; j < result->nrows ; j++ )
	{
	for ( i = 0 ; i < result->ncols ; i++ )
		{
		*resbase = *mat1base + *mat2base ;
		resbase++ ; mat1base++ ; mat2base++ ;
		/* *(resbase++) = *(mat1base++) + *(mat2base++ ); */
		}
	}
if ( is_ephemeral( mat1 ) ) destroy_matrix( mat1 );
if ( is_ephemeral( mat2 ) ) destroy_matrix( mat2 );
return result;
}

MATRIX *matsub( mat1, mat2 )
MATRIX *mat1, *mat2;
{
MATRIX *result;
double *mat1base, *mat2base, *resbase;
int i, j, nlen;
if ( ( mat1->ncols != mat2->ncols ) || ( mat1->nrows != mat2->nrows ) )
	{
	fprintf(stderr,"matsub: args (%dx%d) + (%dx%d) don't conform.\n",
		mat1->nrows, mat1->ncols, mat2->nrows,
		mat2->ncols );
	fprintf(stderr,"matsub: fatal error.  exits. \n");
	Seterr_and_terminate(MATSUB_NONCONFORMITY);
	}
result = create_matrix( mat1->nrows , mat1->ncols, EPHEMERAL);
resbase = result->data;
mat1base = mat1->data;
mat2base = mat2->data;
for ( j = 0 ; j < result->nrows ; j++ )
	{
	for ( i = 0 ; i < result->ncols ; i++ )
		{
		*resbase = *mat1base - *mat2base ;
		resbase++ ; mat1base++ ; mat2base++ ;
		/* *(resbase++) = *(mat1base++) - *(mat2base++ ); */
		}
	}
if ( is_ephemeral( mat1 ) ) destroy_matrix( mat1 );
if ( is_ephemeral( mat2 ) ) destroy_matrix( mat2 );
return result;
}

MATRIX *matmult( mat1, mat2 )
MATRIX *mat1, *mat2;
{
double *mat1base, *mat1loc, *mat2base, *mat2loc, *resbase;
MATRIX *result;
int i, rows, j, nlen;

if ( mat1->ncols != mat2->nrows )
	{
	fprintf(stderr,"matmult: args (%dx%d) * (%dx%d) don't conform.\n",
		mat1->nrows, mat1->ncols, mat2->nrows,
		mat2->ncols );
	fprintf(stderr,"matmult: fatal error.  exits. \n");
	Seterr_and_terminate(MATMULT_NONCONFORMITY);
	}

result = create_matrix( mat1->nrows , mat2->ncols , EPHEMERAL);

resbase = result->data;
mat1base = mat1->data;
mat2base = mat2->data;

for ( j = 0 ; j < result->nrows ; j++ )
	{
	for ( i = 0 ; i < result->ncols ; i++ )
		{
		mat1loc = mat1base;
		mat2loc = mat2base;
		for ( rows = 0 ; rows < mat2->nrows ; rows++ )
			{
			*resbase += *(mat1loc++) * *mat2loc;
			mat2loc += mat2->ncols;
			}
		++resbase;
		++mat2base;
		}
	mat1base += mat1->ncols;
	mat2base = mat2->data;
	}
if ( is_ephemeral( mat1 ) ) destroy_matrix( mat1 );
if ( is_ephemeral( mat2 ) ) destroy_matrix( mat2 );
return result;
}
 
MATRIX *matxdiagasvec( x, d ) MATRIX *x, *d;
{
MATRIX *ans;
double tmp;
int xr, xc, i, j, dr, dc;
xr = x->nrows;
xc = x->ncols;
dr = d->nrows;
dc = d->ncols;

ans = create_matrix( xr, xc, EPHEMERAL );

if ( dc != 1 ) fprintf( stderr, "M+-: matxdiagasvec: d is not a vec.\n");
if ( xc != dr ) fprintf( stderr, "M+-: matxdiagasvec: x and d do not conform\n");

for ( i = 0 ; i < xr ; i++ )
	{
	for ( j = 0 ; j < xc ; j++ )
		{
		MEL( ans, i, j ) = MEL( x, i, j ) * MEL( d, j, 0 );
		}
	}

/* added cleaning 94/11/09 -- pjh */
free_if_ephemeral( x );
free_if_ephemeral( d );

return ans;
}


MATRIX *px1_times_pxq( px1, pxq) /* mult elements of a colvec */
				/* across corresp row of mat */
MATRIX *px1, *pxq;
{
MATRIX *tmp;
double *load, colel;
int i, j;

if ( px1->ncols != 1 )
	{
	fprintf( stderr,"M+-: px1_times_pxq: arg1 not a col-vec. Dies.\n");
	Seterr_and_terminate(PX1XPXQ_ARG1_BAD);
	}
if ( px1->nrows != pxq->nrows )
	{
	fprintf( stderr,"M+-: px1_times_pxq: args not conforming.  Dies.\n");
	Seterr_and_terminate(PX1XPXQ_CONFORMITY);
	}
tmp = matcopy( pxq );
load = tmp->data;
for ( i = 0 ; i < tmp->nrows ; i++ )
	{
	colel = MEL( px1, i, 0);
	for ( j = 0 ; j < tmp->ncols ; j++ )
		{
		*load *= colel ;
		load++ ;
		}
	}
free_if_ephemeral(px1);
free_if_ephemeral(pxq);
return tmp;
}

MATRIX *pxq_divby_px1( pxq, px1) /* divide elements of a colvec */
				/* into corresp row of mat */
MATRIX *px1, *pxq;
{
MATRIX *tmp;
double *load, colel;
int i, j;
if ( px1->ncols != 1 )
	{
	fprintf( stderr,"M+-: pxq_divby_px1: arg2 not a col-vec. Dies.\n");
	Seterr_and_terminate(PXQDPX1_ARG1_BAD);
	}
if ( px1->nrows != pxq->nrows )
	{
	fprintf( stderr,"M+-: pxq_divby_px1: args not conforming.  Dies.\n");
	Seterr_and_terminate(PXQDPX1_CONFORMITY);
	}

tmp = matcopy( pxq );
load = tmp->data;
for (  i = 0 ; i < tmp->nrows ; i++ )
	{
	colel = MEL( px1, i, 0);
	for ( j = 0 ; j < tmp->ncols ; j++ )
		{
		*load = (*load) / colel ;
		load++ ;
		}
	}
free_if_ephemeral(px1);
free_if_ephemeral(pxq);
return tmp;
}
 
MATRIX *scalar_times_matrix( a , X )
double a;
MATRIX *X;
{
MATRIX *tmp;
double *tbase;
int i, nelem;
tmp = matcopy(X);
nelem = get_nelem(tmp);
tbase = tmp->data;
for ( i = 0 ; i < nelem ; i++ ) {
	*tbase *= a ;
	tbase++ ;
}
free_if_ephemeral( X );
return tmp;
}

 
 


void matdump(mat) /* Changed from *matdump to matdump on 4/11/2007: Bryan Comstock */
MATRIX *mat;
{
double *curel;
int outtok = 0;
int nel;

nel = mat->nrows * mat->ncols;

for ( curel = mat->data ;  curel < mat->data + nel ; curel++ )
	{
	printf(  ((fabs(*curel)<.00001) && (fabs(*curel)>0.)) ? "%.4le%c" : "%.4lf%c" , *curel, 
			( outtok++%mat->ncols == mat->ncols-1 )
			? '\n' : ' ' );
	}
/* DOES NOT CLEAN */
}

void fmatdump( of, mat )
MATRIX *mat;
FILE *of;
{
double *curel;
int outtok = 0;
int nel;

nel = mat->nrows * mat->ncols;

for ( curel = mat->data ;  curel < mat->data + nel ; curel++ )
	{
	fprintf( of, ((fabs(*curel)<.00001) && (fabs(*curel)>0.)) ? "%.4le%c" : "%.4lf%c" , *curel,
			( outtok++%mat->ncols == mat->ncols-1 )
			? '\n' : ' ' );
	}
/* DOES NOT CLEAN */
}

void fmatdump_SN( of, mat )
MATRIX *mat;
FILE *of;
{
double *curel;
int outtok = 0;
int nel;

nel = mat->nrows * mat->ncols;

for ( curel = mat->data ;  curel < mat->data + nel ; curel++ )
	{
	fprintf( of,  "%.6le%c" , *curel,
			( outtok++%mat->ncols == mat->ncols-1 )
			? '\n' : ' ' );
	}
/* DOES NOT CLEAN */
}

MATRIX * matread(filename)
char *filename;

{

	MATRIX *X;
	char c, *fmt; 
	FILE *fd, *fopen();
	int wc(), tok, nlines=0, ncols=0, i, j;

	fd = fopen(filename,"r");
	if (fd == NULL)
	{
		printf("\nmatread: couldn't open %s. Retry...\n",filename);
		Seterr_and_terminate(MATREAD_OPEN_FAIL );
	}

	while ((c = fgetc(fd)) && !feof(fd))
	{
		if (c == '\n') nlines++;  /* number of lines */
	}

	fclose(fd);  /* rewind */

	fd = fopen(filename,"r");

	tok = wc(fd);

	fclose(fd);

	if (tok % nlines)  /* test for rectangularity */
	{
		fprintf(stderr,"Unequal #s of numbers per line.\n");
		Seterr_and_terminate(MATREAD_NOT_RECTANGLE);
	}

	ncols = tok/nlines;

	X = create_matrix(nlines,ncols,PERMANENT);

	fd = fopen(filename,"r");

	for ( i = 0 ; i < X->nrows ; i++ )
		{
		for ( j = 0 ; j < X->ncols ; j++ )
			{
			fscanf(fd,"%lf",ELREF(X,i,j));
			}
		}
	return X;
}

 
#define from_S( Sdblptr , Srowintptr , Scolintptr , Matptr ) \
Matptr = create_matrix( *Srowintptr, *Scolintptr , EPHEMERAL ); \
{ \
int i, j, Scol, Srow; \
double *Sload; \
Scol = *Scolintptr; \
Srow = *Srowintptr; \
Sload = Sdblptr; \
for ( j = 0 ; j < Scol ; j++ ) \
	{ \
	for ( i = 0 ; i < Srow ; i++ ) \
		{ \
		MEL( Matptr , i , j ) = (double) * ( Sload ++ ); \
		} \
	} \
}
/* end define |from_S| */

#define to_S( Matptr, Sdblptr ) \
{ \
int i, j; \
double *Sload; \
Sload = Sdblptr; \
for ( j = 0 ; j < Matptr->ncols ; j++ ) \
	{ \
	for ( i = 0 ; i < Matptr->nrows ; i++ ) \
		{ \
		* ( Sload ++ ) = MEL( Matptr , i , j ); \
		} \
	} \
}
 
MATRIX *luinv(X)
MATRIX *X;
{
  MATRIX *Y;
  doublereal det[2] , *work;
  integer job[1];
  int dgefaXXY_(), dgediXXY_(), i;

  doublereal *y;
  integer outsub1, outsub2;
  integer info[1], *ipvt;
  integer nrows, ncols;

  det[0] = 0.;
  det[1] = 0.;

  info[0] = 0;
  job[0] = 0;

  Y = matcopy(X);  /* inversion in situ */

  nrows = X->nrows;
  ncols = X->ncols;

  ipvt = (integer *)malloc((unsigned)nrows*sizeof(integer));
  work = (doublereal *)malloc((unsigned)nrows*sizeof(doublereal));

  y = (doublereal *)Y->data;
  outsub1 = (integer)dgefaXXY_(y,&nrows,&ncols,ipvt,info);

  job[0] = 11;
  outsub2 = (integer)dgediXXY_(y,&nrows,&ncols,ipvt,det,work,job);

  free(ipvt);
  free(work);
  free_if_ephemeral(X);

  return Y;
}

MATRIX *covlag( inmat, lag, demean )
MATRIX *inmat;
int lag, demean;
{
MATRIX *xrows[MAX_COVLAG], *res, *temp;
int n, i, j, nv, q;
double nrec;

n = inmat->nrows;
nrec = (double)1./(double)n;
if ( n > MAX_COVLAG )
	{
	fprintf(stderr,"covlag: arg has > MAX_COVLAG rows. Dies.\n");
	Seterr_and_terminate(EXCEED_MAX_COVLAG);
	}

nv = inmat->ncols;

res = create_matrix( nv, lag*nv, EPHEMERAL );

for ( q = 0 ; q < n ; q++ )
	{
	xrows[q] = extract_rows( inmat, q , q );
	make_permanent(xrows[q]);
	}


for ( i = 0 ; i < lag ; i++ )
	{
	temp = create_matrix( nv, nv, EPHEMERAL );
	for ( j = i ; j < n ; j++ )
		{
		if ( (j-i) < n ) temp = matadd( temp,
			matmult(transp(xrows[j]),xrows[j-i]));
		}
	plug( scalar_times_matrix( nrec, temp ), res, 0, i*nv );
	}

for ( q = 0 ; q < n ; q++ )
	{
	destroy_matrix(xrows[q]);
	}
return res;
}


MATRIX * sweep(mat)

/* algorithm follows Goodnight, Amer Stat Aug 1979 p.149 */

MATRIX * mat;

{

int k , j , i ;
double d , b, *in, *out;
MATRIX *temp, *create_matrix();

temp = create_matrix(mat->nrows,mat->ncols,EPHEMERAL);

in = mat->data;
out = temp->data;

for ( i = 0 ; i < mat->nrows*mat->ncols ; /* copy to out */
		*(out++) = *(in++), i++ )
		;


for ( k = 0 ; k < mat->nrows ; k++ )   /* do 3 */

	{

	d = *(ELREF(temp,k,k));

	for ( j = 0 ; j < mat->nrows ; j++ )   /*do 1 */

		{

		*(ELREF(temp,k,j)) = *(ELREF(temp,k,j)) / d ;

		}   /* 1 */

	for ( i = 0 ; i < mat->nrows ; i++ )  /* do 2 */

		{

		if ( i == k ) goto blast;

		else

			{

			b = *(ELREF(temp,i,k));

			for ( j = 0 ; j < mat->nrows ; j++ )

				{

				*(ELREF(temp,i,j)) =
					*(ELREF(temp,i,j)) - b* *(ELREF(temp,k,j));

				}

			*(ELREF(temp,i,k)) = -b / d;

			}   /* end else: 2 */  

blast: ;		
		}					/* 2 */	

		*(ELREF(temp,k,k)) = (double) 1.0 / d ;

	}  /* 3 */

free_if_ephemeral(mat);
return temp;

}


MATRIX *ident( ord )
int ord;
{
MATRIX *I;
int i;

I = create_matrix( ord, ord, EPHEMERAL );
for ( i = 0 ; i < ord ; i++ )
	*(ELREF(I,i,i)) = (double)1.0;
return I;
}

MATRIX *col_1s( k )
int k;
{
MATRIX *tmp;
int i;
tmp = create_matrix( k , 1 , EPHEMERAL );
for ( i = 0 ; i < k ; i++ )
	{
	MEL(tmp,i,0) = 1.;
	}
return tmp;
}

int wc(file)

FILE *file;

{
	int c,nl,nw,nc,inword;

	inword = 0;
	nl = nw = nc = 0;
	while ((c = fgetc(file)) && !feof(file))
	{
		++nc;
		if (c == '\n') ++nl;
		if (c == ' ' || c == '\n' || c == '\t')
			inword = 0;
		else if (inword == 0)
		{
			inword = 1;
			++nw;
		}
	}
	return nw;
}

int nchanges(X)
MATRIX *X;
{
/* returns integer telling how often the value of X */
/* changes from row to row.  X must be column vector */

int tmp = 1, istart, i;

if (X->ncols != 1)
	{
	fprintf(stderr,"nchanges:  must be column vector; ncols = %d.\n",
				X->ncols);
	fprintf(stderr,"nchanges: exiting.\n");
	exit(1);
	}

istart = MEL( X , 0 , 0 );

for ( i = 1 ; i < X->nrows ; i++ )
	{
	if ( MEL ( X , i , 0 ) != istart )
		{
		tmp++;
		istart = MEL ( X , i , 0 );
		}
	}
return tmp;
}
 


MATRIX *Cchol( x )      /* returns Lt (transpose of L -- lower triangular */
MATRIX *x;              /* where   L Lt = x                               */
{                       /*                   L = x 0 0 0                  */
int p, i, j, k;         /*                       x x 0 0                  */
MATRIX *tmp;            /*                       x x x 0                  */
double *tmphead;        /*                       x x x x                  */

p = x->nrows;

if ( p != x->ncols )
  {
        fprintf(stderr,"Cchol: not a square matrix\n. Dies.\n");
        Seterr_and_terminate(CCHOL_NOT_SQUARE);
      }

tmp = matcopy(x);
tmphead = tmp->data;

#define loctmp( i , j ) *(tmphead + ( i * p ) + j )

for ( i = 0 ; i < p ; i++ )
  {
        double accum = 0.;

        for ( k = 0 ; k <= i-1 ; k++ )
	  {
                accum += MEL(tmp,k,i) * MEL(tmp,k,i);
	      }

        if ( ( MEL(tmp,i,i) - accum ) >= 0. )
	  {
                loctmp( i , i ) = (double) sqrt( MEL(tmp,i,i) - accum );
	      }
        else
	  {
                fprintf( stderr , "Cchol: illegal diag element encountered.");
                fprintf( stderr , "Dies.\n");
                Seterr_and_terminate(CCHOL_FAIL);
	      }


                                                      
        for ( j = i+1  ; j < p ; j++ )
	  {
                double accum2 = 0.;
                for ( k = 0 ; k <= i-1 ; k++ )
		  {
                        accum2 += MEL(tmp,k,i) * MEL(tmp,k,j);
		      }
                if ( MEL(tmp,i,i) > 0. ) loctmp( i , j ) = ( MEL(tmp,i,j) - 
						      accum2 )/MEL(tmp,i,i);
                else
		  {
                        fprintf( stderr,"Cchol: non-positive pivot.  Dies.\n");
                        Seterr_and_terminate(CCHOL_FAIL);
		      }
                loctmp( j , i ) = 0.;
	      }
      }
if ( is_ephemeral(x) ) destroy_matrix(x) ;

return tmp;

}

/* This point marks the end of file chanmat.c :: B.A.C. */
/********************************************************/
