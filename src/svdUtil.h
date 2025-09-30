#ifndef RF_SVD_UTIL_H
#define RF_SVD_UTIL_H
void svdcmp(double **a, int m, int n, double ***uptr, double **wptr, double ***vptr);
char svdchk(double **a, uint m, uint n, double **u, double *w, double **v);
double **svdinv(double **u, double *w, double **v, uint m, uint n, uint singularity);
void free_svdcmp(double **a, int m, int n, double **u, double *w, double **v);
void svbksb(double **u, double *w, double **v, uint m, uint n, double *b, double *x);
double **matrixCopy(double **a, uint m, uint n);
double **matrixTrans(double **a, uint m, uint n);
double **matrixMult(double **a, double **b, uint m, uint n, uint p);
void matrixPrint(double **x, uint m, uint n);
double pythag(double a, double b);
void harness(void);
#endif
