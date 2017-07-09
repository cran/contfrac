/* This generated by running tools::package_native_routine_registration_skeleton(".") */

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void c_contfrac(void *, void *, void *, void *, void *);
extern void c_contfrac_complex(void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_convergents(void *, void *, void *, void *, void *, void *);
extern void c_convergents_complex(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"c_contfrac",            (DL_FUNC) &c_contfrac,             5},
    {"c_contfrac_complex",    (DL_FUNC) &c_contfrac_complex,     8},
    {"c_convergents",         (DL_FUNC) &c_convergents,          6},
    {"c_convergents_complex", (DL_FUNC) &c_convergents_complex, 11},
    {NULL, NULL, 0}
};

void R_init_contfrac(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}