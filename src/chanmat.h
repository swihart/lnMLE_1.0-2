/* gee support @(#) chanmat.h 3.1 94/03/03 */
#include "chanmatstruct.h"
#include "chanmatfuns.h"
#include <setjmp.h>
/* extern jmp_buf env;

/* Bryan commented this out... it exists in chanmatstruct.h and is wreaking havoc on my program */
/* #define errorbranch( exitlabel ) \ */
/* if ( setjmp(env) != 0 ) \ */
/*        { \ */
/*        fprintf(stderr,"chanmat error detected, returning to caller\n"); */
/*        goto exitlabel ; */
/*        } */

