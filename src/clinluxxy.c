
/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
/* From this point onward is the file clinluxxy.c :: B.A.C. */
#include "f2c.h"

doublereal dasumXXY_(n, dx, incx)
integer *n;
doublereal *dx;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i, m, ns, mp1;

/* ***BEGIN PROLOGUE  DASUM */
/* ***DATE WRITTEN   791001   (YYMMDD) */
/* ***REVISION DATE  820801   (YYMMDD) */
/* ***CATEGORY NO.  D1A3A */
/* ***KEYWORDS  ADD,BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAGNITUDE,SUM, */
/*             VECTOR */
/* ***AUTHOR  LAWSON, C. L., (JPL) */
/*           HANSON, R. J., (SNLA) */
/*           KINCAID, D. R., (U. OF TEXAS) */
/*           KROGH, F. T., (JPL) */
/* ***PURPOSE  Sum of magnitudes of d.p. vector components */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */

/*     --Output-- */
/*    DASUM  double precision result (zero if N .LE. 0) */

/*     Returns sum of magnitudes of double precision DX. */
/*     DASUM = sum from 0 to N-1 of DABS(DX(1+I*INCX)) */
/* ***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, 
*/
/*                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL */
/*                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 
*/
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  DASUM */

/* ***FIRST EXECUTABLE STATEMENT  DASUM */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENTS NOT EQUAL TO 1. */

    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
	if (dx[i] < 0.) {
	    dx[i] = -dx[i];
	}
/*          DASUM = DASUM + DABS(DX(I)) */
	ret_val += dx[i];
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENTS EQUAL TO 1. */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6. */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= i__2; ++i) {
	if (dx[i] < 0.) {
	    dx[i] = -dx[i];
	}
/*         DASUM = DASUM + DABS(DX(I)) */
	ret_val += dx[i];
/* L30: */
    }
    if (*n < 6) {
	return ret_val;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= i__2; i += 6) {
	ret_val = ret_val + (d__1 = dx[i], abs(d__1)) + (d__2 = dx[i + 1], 
		abs(d__2)) + (d__3 = dx[i + 2], abs(d__3)) + (d__4 = dx[i + 3]
		, abs(d__4)) + (d__5 = dx[i + 4], abs(d__5)) + (d__6 = dx[i + 
		5], abs(d__6));
/* L50: */
    }
    return ret_val;
} /* dasum_ */

/* daxpy.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Subroutine */ int daxpyXXY_(n, da, dx, incx, dy, incy)
integer *n;
doublereal *da, *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dy[i] += *da * dx[i];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 4) {
	dy[i] += *da * dx[i];
	dy[i + 1] += *da * dx[i + 1];
	dy[i + 2] += *da * dx[i + 2];
	dy[i + 3] += *da * dx[i + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */

/* ddot.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


doublereal ddotXXY_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static integer i, m, ix, iy, ns, mp1;

/* ***BEGIN PROLOGUE  DDOT */
/* ***DATE WRITTEN   791001   (YYMMDD) */
/* ***REVISION DATE  820801   (YYMMDD) */
/* ***CATEGORY NO.  D1A4 */
/* ***KEYWORDS  BLAS,DOUBLE PRECISION,INNER PRODUCT,LINEAR ALGEBRA,VECTOR 
*/
/* ***AUTHOR  LAWSON, C. L., (JPL) */
/*           HANSON, R. J., (SNLA) */
/*           KINCAID, D. R., (U. OF TEXAS) */
/*           KROGH, F. T., (JPL) */
/* ***PURPOSE  D.P. inner product of d.p. vectors */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */
/*       DY  double precision vector with N elements */
/*     INCY  storage spacing between elements of DY */

/*     --Output-- */
/*     DDOT  double precision dot product (zero if N .LE. 0) */

/*     Returns the dot product of double precision DX and DY. */
/*     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY) */
/*     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is */
/*     defined in a similar way using INCY. */
/* ***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, 
*/
/*                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL */
/*                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 
*/
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  DDOT */

/* ***FIRST EXECUTABLE STATEMENT  DDOT */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == *incy) {
	if ((i__1 = *incx - 1) < 0) {
	    goto L5;
	} else if (i__1 == 0) {
	    goto L20;
	} else {
	    goto L60;
	}
    }
L5:

/*         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS. */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	ret_val += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return ret_val;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1. */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5. */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	ret_val += dx[i] * dy[i];
/* L30: */
    }
    if (*n < 5) {
	return ret_val;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 5) {
	ret_val = ret_val + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[i + 2] 
		* dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 4] * dy[i + 4];
/* L50: */
    }
    return ret_val;

/*         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
	ret_val += dx[i] * dy[i];
/* L70: */
    }
    return ret_val;
} /* ddot_ */

/* dgeco.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dgecoXXY_(a, lda, n, ipvt, rcond, z)
doublereal *a;
integer *lda, *n, *ipvt;
doublereal *rcond, *z;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign();

    /* Local variables */
    extern doublereal ddotXXY_();
    static integer info;
    extern /* Subroutine */ int dgefaXXY_();
    static integer j, k, l;
    static doublereal s, t;
    extern /* Subroutine */ int dscalXXY_();
    extern doublereal dasumXXY_();
    static doublereal anorm;
    extern /* Subroutine */ int daxpyXXY_();
    static doublereal ynorm;
    static integer kb;
    static doublereal ek, sm, wk;
    static integer kp1;
    static doublereal wkm;

/* ***BEGIN PROLOGUE  DGECO */
/* ***DATE WRITTEN   780814   (YYMMDD) */
/* ***REVISION DATE  820801   (YYMMDD) */
/* ***CATEGORY NO.  D2A1 */
/* ***KEYWORDS  CONDITION,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK, 
*/
/*             MATRIX */
/* ***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO) */
/* ***PURPOSE  Factors a double precision matrix by Gaussian elimination 
*/
/*            and estimates the condition of the matrix. */
/* ***DESCRIPTION */

/*     DGECO factors a double precision matrix by Gaussian elimination */
/*     and estimates the condition of the matrix. */

/*     If  RCOND  is not needed, DGEFA is slightly faster. */
/*     To solve  A*X = B , follow DGECO by DGESL. */
/*     To compute  INVERSE(A)*C , follow DGECO by DGESL. */
/*     To compute  DETERMINANT(A) , follow DGECO by DGEDI. */
/*     To compute  INVERSE(A) , follow DGECO by DGEDI. */

/*     On Entry */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                the matrix to be factored. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        A       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                The factorization can be written  A = L*U  where */
/*                L  is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        IPVT    INTEGER(N) */
/*                an INTEGER vector of pivot indices. */

/*        RCOND   DOUBLE PRECISION */
/*                an estimate of the reciprocal condition of  A . */
/*                For the system  A*X = B , relative perturbations */
/*                in  A  and  B  of size  EPSILON  may cause */
/*                relative perturbations in  X  of size  EPSILON/RCOND . 
*/
/*                If  RCOND  is so small that the logical expression */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                is true, then  A  may be singular to working */
/*                precision.  In particular,  RCOND  is zero  if */
/*                exact singularity is detected or the estimate */
/*                underflows. */

/*        Z       DOUBLE PRECISION(N) */
/*                a work vector whose contents are usually unimportant. */
/*                If  A  is close to a singular matrix, then  Z  is */
/*                an approximate null vector in the sense that */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     LINPACK.  This version dated 08/14/78 . */
/*     Cleve Moler, University of New Mexico, Argonne National Lab. */

/*     Subroutines and Functions */

/*     LINPACK DGEFA */
/*     BLAS DAXPY,DDOT,DSCAL,DASUM */
/*     Fortran DABS,DMAX1,DSIGN */
/* ***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W., */
/*                 *LINPACK USERS  GUIDE*, SIAM, 1979. */
/* ***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGEFA,DSCAL */
/* ***END PROLOGUE  DGECO */


/*     COMPUTE 1-NORM OF A */

/* ***FIRST EXECUTABLE STATEMENT  DGECO */
    /* Parameter adjustments */
    --z;
    --ipvt;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = dasumXXY_(n, &a[j * a_dim1 + 1], &c__1);
	anorm = max(d__1,d__2);
/* L10: */
    }

/*     FACTOR */

    dgefaXXY_(&a[a_offset], lda, n, &ipvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E . */
/*     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE */
/*     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE */
/*     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID */
/*     OVERFLOW. */

/*     SOLVE TRANS(U)*W = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z[j] = 0.;
/* L20: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z[k] != 0.) {
	    d__1 = -z[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z[k], abs(d__1)) <= (d__2 = a[k + k * a_dim1], abs(
		d__2))) {
	    goto L30;
	}
	s = (d__1 = a[k + k * a_dim1], abs(d__1)) / (d__2 = ek - z[k], abs(
		d__2));
	dscalXXY_(n, &s, &z[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z[k];
	wkm = -ek - z[k];
	s = abs(wk);
	sm = abs(wkm);
	if (a[k + k * a_dim1] == 0.) {
	    goto L40;
	}
	wk /= a[k + k * a_dim1];
	wkm /= a[k + k * a_dim1];
	goto L50;
L40:
	wk = 1.;
	wkm = 1.;
L50:
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    sm += (d__1 = z[j] + wkm * a[k + j * a_dim1], abs(d__1));
	    z[j] += wk * a[k + j * a_dim1];
	    s += (d__1 = z[j], abs(d__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	t = wkm - wk;
	wk = wkm;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    z[j] += t * a[k + j * a_dim1];
/* L70: */
	}
L80:
L90:
	z[k] = wk;
/* L100: */
    }
    s = 1. / dasumXXY_(n, &z[1], &c__1);
    dscalXXY_(n, &s, &z[1], &c__1);

/*     SOLVE TRANS(L)*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if (k < *n) {
	    i__2 = *n - k;
	    z[k] += ddotXXY_(&i__2, &a[k + 1 + k * a_dim1], &c__1, &z[k + 1], &
		    c__1);
	}
	if ((d__1 = z[k], abs(d__1)) <= 1.) {
	    goto L110;
	}
	s = 1. / (d__1 = z[k], abs(d__1));
	dscalXXY_(n, &s, &z[1], &c__1);
L110:
	l = ipvt[k];
	t = z[l];
	z[l] = z[k];
	z[k] = t;
/* L120: */
    }
    s = 1. / dasumXXY_(n, &z[1], &c__1);
    dscalXXY_(n, &s, &z[1], &c__1);

    ynorm = 1.;

/*     SOLVE L*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = z[l];
	z[l] = z[k];
	z[k] = t;
	if (k < *n) {
	    i__2 = *n - k;
	    daxpyXXY_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &z[k + 1], &c__1)
		    ;
	}
	if ((d__1 = z[k], abs(d__1)) <= 1.) {
	    goto L130;
	}
	s = 1. / (d__1 = z[k], abs(d__1));
	dscalXXY_(n, &s, &z[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1. / dasumXXY_(n, &z[1], &c__1);
    dscalXXY_(n, &s, &z[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE  U*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z[k], abs(d__1)) <= (d__2 = a[k + k * a_dim1], abs(d__2)))
		 {
	    goto L150;
	}
	s = (d__1 = a[k + k * a_dim1], abs(d__1)) / (d__2 = z[k], abs(d__2));
	dscalXXY_(n, &s, &z[1], &c__1);
	ynorm = s * ynorm;
L150:
	if (a[k + k * a_dim1] != 0.) {
	    z[k] /= a[k + k * a_dim1];
	}
	if (a[k + k * a_dim1] == 0.) {
	    z[k] = 1.;
	}
	t = -z[k];
	i__2 = k - 1;
	daxpyXXY_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z[1], &c__1);
/* L160: */
    }
/*     MAKE ZNORM = 1.0 */
    s = 1. / dasumXXY_(n, &z[1], &c__1);
    dscalXXY_(n, &s, &z[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* dgeco_ */

/* dgedi.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */


/* Subroutine */ int dgediXXY_(a, lda, n, ipvt, det, work, job)
doublereal *a;
integer *lda, *n, *ipvt;
doublereal *det, *work;
integer *job;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, k, l;
    static doublereal t;
    extern /* Subroutine */ int dscalXXY_(), dswapXXY_(), daxpyXXY_();
    static integer kb, kp1, nm1;
    static doublereal ten;


/*     dgedi computes the determinant and inverse of a matrix */
/*     using the factors computed by dgeco or dgefa. */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the output from dgeco or dgefa. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        ipvt    integer(n) */
/*                the pivot vector from dgeco or dgefa. */

/*        work    double precision(n) */
/*                work vector.  contents destroyed. */

/*        job     integer */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     on return */

/*        a       inverse of original matrix if requested. */
/*                otherwise unchanged. */

/*        det     double precision(2) */
/*                determinant of original matrix if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. dabs(det(1)) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        it will not occur if the subroutines are called correctly */
/*        and if dgeco has set rcond .gt. 0.0 or dgefa has set */
/*        info .eq. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal,dswap */
/*     fortran dabs,mod */

/*     internal variables */



/*     compute determinant */

    /* Parameter adjustments */
    --work;
    --det;
    --ipvt;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (ipvt[i] != i) {
	    det[1] = -det[1];
	}
	det[1] = a[i + i * a_dim1] * det[1];
/*        ...exit */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (abs(det[1]) >= 1.) {
	    goto L20;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (abs(det[1]) < ten) {
	    goto L40;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     compute inverse(u) */

    if (*job % 10 == 0) {
	goto L150;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
	t = -a[k + k * a_dim1];
	i__2 = k - 1;
	dscalXXY_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[k + j * a_dim1];
	    a[k + j * a_dim1] = 0.;
	    daxpyXXY_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        form inverse(u)*inverse(l) */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L140;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	kp1 = k + 1;
	i__2 = *n;
	for (i = kp1; i <= i__2; ++i) {
	    work[i] = a[i + k * a_dim1];
	    a[i + k * a_dim1] = 0.;
/* L110: */
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = work[j];
	    daxpyXXY_(n, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L120: */
	}
	l = ipvt[k];
	if (l != k) {
	    dswapXXY_(n, &a[k * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1], &c__1);
	}
/* L130: */
    }
L140:
L150:
    return 0;
} /* dgedi_ */

/* dgefa.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */


/* Subroutine */ int dgefaXXY_(a, lda, n, ipvt, info)
doublereal *a;
integer *lda, *n, *ipvt, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, l;
    static doublereal t;
    extern /* Subroutine */ int dscalXXY_(), daxpyXXY_();
    extern integer idamaxXXY_();
    static integer kp1, nm1;

/*     dgefa factors a double precision matrix by gaussian elimination. */

/*     dgefa is usually called by dgeco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for dgeco) = (1 + 9/n)*(time for dgefa) . */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the matrix to be factored. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    integer(n) */
/*                an integer vector of pivot indices. */

/*        info    integer */
/*                = 0  normal value. */
/*                = k  if  u(k,k) .eq. 0.0 .  this is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that dgesl or dgedi will divide by zero */
/*                     if called.  use  rcond  in dgeco for a reliable */
/*                     indication of singularity. */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal,idamax */

/*     internal variables */



/*     gaussian elimination with partial pivoting */

    /* Parameter adjustments */
    --ipvt;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        find l = pivot index */

	i__2 = *n - k + 1;
	l = idamaxXXY_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        zero pivot implies this column already triangularized */

	if (a[l + k * a_dim1] == 0.) {
	    goto L40;
	}

/*           interchange if necessary */

	if (l == k) {
	    goto L10;
	}
	t = a[l + k * a_dim1];
	a[l + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L10:

/*           compute multipliers */

	t = -1. / a[k + k * a_dim1];
	i__2 = *n - k;
	dscalXXY_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           row elimination with column indexing */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[l + j * a_dim1];
	    if (l == k) {
		goto L20;
	    }
	    a[l + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
L20:
	    i__3 = *n - k;
	    daxpyXXY_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * 
		    a_dim1], &c__1);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    if (a[*n + *n * a_dim1] == 0.) {
	*info = *n;
    }
    return 0;
} /* dgefa_ */

/* dgidamax.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


integer idamaxXXY_(n, dx, incx)
integer *n;
doublereal *dx;
integer *incx;
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal dmax_;
    static integer i, ix;


/*     finds the index of element having max. absolute value. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    dmax_ = abs(dx[1]);
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if ((d__1 = dx[ix], abs(d__1)) <= dmax_) {
	    goto L5;
	}
	ret_val = i;
	dmax_ = (d__1 = dx[ix], abs(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    dmax_ = abs(dx[1]);
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if ((d__1 = dx[i], abs(d__1)) <= dmax_) {
	    goto L30;
	}
	ret_val = i;
	dmax_ = (d__1 = dx[i], abs(d__1));
L30:
	;
    }
    return ret_val;
} /* idamax_ */

/* dscal.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Subroutine */ int dscalXXY_(n, da, dx, incx)
integer *n;
doublereal *da, *dx;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*     scales a vector by a constant. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
	dx[i] = *da * dx[i];
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= i__2; ++i) {
	dx[i] = *da * dx[i];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= i__2; i += 5) {
	dx[i] = *da * dx[i];
	dx[i + 1] = *da * dx[i + 1];
	dx[i + 2] = *da * dx[i + 2];
	dx[i + 3] = *da * dx[i + 3];
	dx[i + 4] = *da * dx[i + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

/* dswap.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Subroutine */ int dswapXXY_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     interchanges two vectors. */
/*     uses unrolled loops for increments equal one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dtemp = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */


/*       clean-up loop */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dtemp = dx[i];
	dx[i] = dy[i];
	dy[i] = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 3) {
	dtemp = dx[i];
	dx[i] = dy[i];
	dy[i] = dtemp;
	dtemp = dx[i + 1];
	dx[i + 1] = dy[i + 1];
	dy[i + 1] = dtemp;
	dtemp = dx[i + 2];
	dx[i + 2] = dy[i + 2];
	dy[i + 2] = dtemp;
/* L50: */
    }
    return 0;
} /* dswap_ */


double d_sign(a,b)
doublereal *a, *b;
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}


