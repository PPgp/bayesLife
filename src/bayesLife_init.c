#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void doDL(void *, void *, void *, void *, void *, void *);
extern void dologdensityTrianglekz(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"doDL",                   (DL_FUNC) &doDL,                    6},
    {"dologdensityTrianglekz", (DL_FUNC) &dologdensityTrianglekz, 14},
    {NULL, NULL, 0}
};

void R_init_bayesLife(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}