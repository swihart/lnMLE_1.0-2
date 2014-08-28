/* chanmatfuns.h */
/* gee support @(#) chanmatfuns.h 3.2 94/03/09 */


extern MATRIX *create_matrix(),
     *matread(),
     *matcopy(),
     *extract_rows(),
     *matadd(),
     *matsub(),
     *matmult(),
     *transp(),
     *Cchol(),
     *sweep(),
     *col_1s(),
     *matabs(),
     *matexp(),
     *px1_times_pxq(),
     *pxq_divby_px1(),
     *scalar_times_matrix(),
     *ident(),
     *form_diag(),
     *diaginv(),
     *diagsqrt(),
     *corner(),
     *covlag(),
     *toeplitz(),
     *band(),
     *matcumnorm(),
     *extract_cols(),
     *matxdiagasvec(),
     /* following two functions added by pj catalano  */
     *matnpdf(),
     *matncdf();

extern double matmax(), elsum();

extern void matdump(), plug(), destroy_matrix(), fmatdump();

extern MATRIX *luinv();


