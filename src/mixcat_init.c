#include <stdlib.h> 
#include <R_ext/Rdynload.h>

/* .C calls */
extern void npmltd(int *y, int *q1, int *N1, int *m1, int *CN, int *npk1, int *np1, double *model,
            double *eps1, double *strvlint, double *strvlreg, double *strvlmp, double *strvlm,
            double *out, int *EBind1, double *outEB, int *link1, int *maxit1,
            double *rslp, int *rslpind1, double *outFitted, double *outProb, double *tol1,
            double *npoind, int *T1, double *outparcvmat);

static const R_CMethodDef CEntries[] = {
    {"npmltd",      (DL_FUNC) &npmltd,      26},
    {NULL, NULL, 0}
};

void R_init_BNSP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

