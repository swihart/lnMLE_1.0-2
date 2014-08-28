/* gee support @(#) chanmatstruct.h 3.5 94/05/09 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <setjmp.h>
/* f2c.h  --  Standard Fortran to C header file */
 
/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
 
        - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */
 
#ifndef F2C_INCLUDE
#define F2C_INCLUDE
 
typedef int integer;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
#endif
 

#define PERMANENT 1
#define EPHEMERAL 0
#define MAX_COVLAG 30
#define NO_ERROR 0
#define UNKNOWN_FAILURE 1
#define NO_MEM_MATSTRUCT 2
#define NO_MEM_MATDATA 3
#define SPLIT_FAIL 4
#define MATMULT_NONCONFORMITY 5
#define MATADD_NONCONFORMITY 6
#define CCHOL_FAIL 7
#define CORNER_FAIL 8
#define EXCEED_MAX_COVLAG 9
#define BAD_TOEPLITZ_ARG 10
#define PLUG_FAIL 11
#define PX1XPXQ_ARG1_BAD 12
#define PX1XPXQ_CONFORMITY 13
#define PXQDPX1_ARG1_BAD 14
#define PXQDPX1_CONFORMITY 15
#define CCHOL_NOT_SQUARE 16
#define MATREAD_OPEN_FAIL 17
#define MATREAD_NOT_RECTANGLE 18

 

typedef struct matrix
		{
		int nrows, ncols;
		double *data;
		int permanence;
		} MATRIX;
/* element reference is handled principally by MEL */
#define ELREF( matp , s1, s2 ) ((matp)->data)+(s2)+((s1)*(matp->ncols))
#define MEL(X ,i, j) (*(ELREF( (X), (i), (j) ) ))
/* $ y = MEL(X,i,j) :\Rightarrow y \in {\cal R} \wedge y = x_{ij} $ */
#define get_nelem( x ) (((x)->nrows) * ((x)->ncols))
 

#include <setjmp.h>

static jmp_buf env;

#define Seterr_and_terminate( Code ) { fprintf(stderr, \
           "chanmat library error Code , returning.\n"); longjmp(env,1); }

#define errorbranch( exitlabel ) \
if ( setjmp(env) != 0 ) \
        { \
        fprintf(stderr,"chanmat error detected, returning to caller\n"); \
        goto exitlabel ; \
        }
 
#ifdef S91
#include "/usr/local/S.91-ATT/include/S.h"
#define calloc(x,y) Calloc(x,y)
#define cfree Free
#define free Free
#else
#define calloc(x,y) malloc((unsigned)((x)*(y)))
#define cfree free
#endif

#define is_permanent( x ) (x)->permanence == PERMANENT
#define is_ephemeral( x ) (x)->permanence == EPHEMERAL
#define make_permanent( x ) (x)->permanence = PERMANENT;
#define make_ephemeral( x ) (x)->permanence = EPHEMERAL;

#define free_if_ephemeral( x ) if (is_ephemeral((x))) destroy_matrix((x))
 
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
 

