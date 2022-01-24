#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void C_call_nmf(void *, void *, void *, void *, void *);
extern void C_call_suitor(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"C_call_nmf",    (DL_FUNC) &C_call_nmf,    5},
    {"C_call_suitor", (DL_FUNC) &C_call_suitor, 9},
    {NULL, NULL, 0}
};

void R_init_SUITOR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
