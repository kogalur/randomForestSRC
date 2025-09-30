
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "nrutil.h"
#include "error.h"
unsigned int upower (unsigned int x, unsigned int n) {
  unsigned int p;
  if ((x >= 2) & (n > (sizeof(unsigned int) * 8) - 1)) {
    nrerror("Overflow in upower(), exponent too large.");
  }
  for (p = 1; n > 0; --n) {
    p = p * x;
  }
  return p;
}
unsigned int upower2 (unsigned int n) {
  unsigned int p;
  if (n > (sizeof(unsigned int) * 8) - 1) {
    nrerror("Overflow in upower2(), exponent too large.");
  }
  p = ((unsigned int) 1) << n;
  return p;
}
unsigned int ulog2 (unsigned int n) {
  unsigned int p;
  p = 0;
  while (n > 1) {
    n = n >> 1;
    p++;
  }
  return p;
}
void hpsort(double *ra, unsigned int n) {
  unsigned int i, ir, j, l;
  double rra;
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    }
    else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      }
      else {
        j = ir+1;
      }
    }
    ra[i] = rra;
  }
}
void hpsortui(unsigned int *ra, unsigned int n) {
  unsigned int i, ir, j, l;
  unsigned int rra;
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    }
    else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      }
      else {
        j = ir+1;
      }
    }
    ra[i] = rra;
  }
}
void hpsorti(int *ra, unsigned int n) {
  unsigned int i, ir, j, l;
  int rra;
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    }
    else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      }
      else {
        j = ir+1;
      }
    }
    ra[i] = rra;
  }
}
#ifdef SWAP
#undef SWAP
#endif
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50
void qksort(double *arr, unsigned int n) {
  unsigned int i, j, k, l;
  unsigned int ir;
  unsigned int *istack, jstack;
  double a, temp;
  if (n < 1) nrerror("\n n of zero (0) length in indexx().");
  l  = 1;
  ir = n;
  jstack = 0;
  istack = uivector(1, NSTACK);
  for (;;) {
    if (ir-l < M) {
      for (j = l+1; j <= ir; j++) {
        a=arr[j];
        for (i = j-1; i >= l; i--) {
          if (arr[i] <= a) break;
          arr[i+1] = arr[i];
        }
        arr[i+1] = a;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--]; 
    } else {
      k = (l+ir) >> 1; 
      SWAP(arr[k], arr[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l], arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1], arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l], arr[l+1]);
      }
      i = l+1; 
      j = ir;
      a = arr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i], arr[j]);
      }
      arr[l+1] = arr[j]; 
      arr[j] = a;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in sort().");
      if (ir-i+1 >= j-l) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j-1;
      }
      else {
        istack[jstack] = j-1;
        istack[jstack-1] = l;
        l=i;
      }
    }
  }
  free_uivector(istack,1,NSTACK);
}
#undef SWAP
#undef M
#undef NSTACK
#ifdef SWAP
#undef SWAP
#endif
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
void indexx(unsigned int n, double *arr, unsigned int *indx) {
  unsigned int i, j, k, l;
  unsigned int indxt, itemp, ir;
  unsigned int *istack, jstack;
  double a;
  if (n < 1) nrerror("\n n of zero (0) length in indexx().");
  l  = 1;
  ir = n;
  jstack = 0;
  istack = uivector(1, NSTACK);
  for (j=1; j<=n; j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j = l+1; j <= ir; j++) {
        indxt = indx[j];
        a = arr[indxt];
        for (i=j-1; i>=l; i--) {
          if (arr[indx[i]] <= a) break;
          indx[i+1] = indx[i];
        }
        indx[i+1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--];
    }
    else {
      k = (l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if (arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l], indx[ir])
      }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP(indx[l+1], indx[ir])
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
        SWAP(indx[l], indx[l+1])
      }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      a = arr[indxt];
      for (;;) {
        do i++; while (arr[indx[i]] < a);
        do j--; while (arr[indx[j]] > a);
        if (j < i) break;
        SWAP(indx[i], indx[j])
      }
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in indexx().");
      if (ir-i+1 >= j-l) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j-1;
      }
      else {
        istack[jstack] = j-1;
        istack[jstack-1] = l;
        l = i;
      }
    }
  }
  free_uivector(istack, 1, NSTACK);
}
#undef SWAP
#undef M
#undef NSTACK
void nrerror(char error_text[]) {
  RF_nativeError("\nRF-SRC");
  RF_nativeError("\nRF-SRC:  *** ERROR *** ");
  RF_nativeError("\nRF-SRC:  Numerical Recipes Run-Time Error:");
  RF_nativeError("\nRF-SRC:  %s", error_text);
  RF_nativeError("\nRF-SRC:  Please Contact Technical Support.");
  RF_nativeExit();
}
void *gblock(size_t size) {
  void *v = (void *) malloc(size);
  if (!v) nrerror("\n  Allocation Failure in gblock().");
  return v;
}
void free_gblock(void *v, size_t size) {
  free((FREE_ARG) v);
}
void *gvector(unsigned long long nl, unsigned long long nh, size_t size) {
  if (nh < nl) nrerror("\n  Illegal indices in gvector().");
  void *v = gblock((size_t) ((nh-nl+1+NR_END) * size));
  return v;
}
void free_gvector(void *v, unsigned long long nl, unsigned long long nh, size_t size) {
  if (nh < nl) nrerror("\n  Illegal indices in free_gvector().");
  free_gblock(v, (nh-nl+1+NR_END) * size);
}
char *cvector(unsigned long long nl, unsigned long long nh) {
  return ((char *) gvector(nl, nh, sizeof(char)) -nl+NR_END);
}
void free_cvector(char *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(char));
}
char **cmatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  char **v = (char **) new_vvector(nrl, nrh, NRUTIL_CPTR);
  for(unsigned long long i = nrl; i <= nrh; i++) {
    v[i] = cvector(ncl, nch);
  }
  return v;
}
void free_cmatrix(char **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = nrl; i <= nrh; i++) {
    free_cvector(v[i], ncl, nch);
  }
  free_new_vvector(v, nrl, nrh, NRUTIL_CPTR);
}
int *ivector(unsigned long long nl, unsigned long long nh) {
  return ((int *) gvector(nl, nh, sizeof(int)) -nl+NR_END);
}
void free_ivector(int *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(int));
}
int **imatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  int **v = (int **) new_vvector(nrl, nrh, NRUTIL_IPTR);
  for(unsigned long long i = nrl; i <= nrh; i++) {
    v[i] = ivector(ncl, nch);
  }
  return v;
}
void free_imatrix(int **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = nrl; i <= nrh; i++) {
    free_ivector(v[i], ncl, nch);
  }
  free_new_vvector(v, nrl, nrh, NRUTIL_IPTR);
}
unsigned int *uivector(unsigned long long nl, unsigned long long nh) {
  return ((unsigned int *) gvector(nl, nh, sizeof(unsigned int)) -nl+NR_END);
}
void free_uivector(unsigned int *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(unsigned int));
}
unsigned int **uimatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  unsigned int **v = (unsigned int **) new_vvector(nrl, nrh, NRUTIL_UPTR);
  for(unsigned long long i = nrl; i <= nrh; i++) {
    v[i] = uivector(ncl, nch);
  }
  return v;
}
void free_uimatrix(unsigned int **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = nrl; i <= nrh; i++) {
    free_uivector(v[i], ncl, nch);
  }
  free_new_vvector(v, nrl, nrh, NRUTIL_UPTR);
}
unsigned long *ulvector(unsigned long long nl, unsigned long long nh) {
  return ((unsigned long *) gvector(nl, nh, sizeof(unsigned long)) -nl+NR_END);
}
void free_ulvector(unsigned long *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(unsigned long));
}
double *dvector(unsigned long long nl, unsigned long long nh) {
  return ((double *) gvector(nl, nh, sizeof(double)) -nl+NR_END);
}
void free_dvector(double *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(double));
}
double **dmatrix(unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  double **v = (double **) new_vvector(nrl, nrh, NRUTIL_DPTR);
  for(unsigned long long i = nrl; i <= nrh; i++) {
    v[i] = dvector(ncl, nch);
  }
  return v;
}
void free_dmatrix(double **v, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = nrl; i <= nrh; i++) {
    free_dvector(v[i], ncl, nch);
  }
  free_new_vvector(v, nrl, nrh, NRUTIL_DPTR);
}
double ***dmatrix3(unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  double ***v = (double ***) new_vvector(n3l, n3h, NRUTIL_DPTR2);
  for(unsigned long long i = n3l; i <= n3h; i++) {
    v[i] = dmatrix(nrl, nrh, ncl, nch);
  }
  return v;
}
void free_dmatrix3(double ***v, unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = n3l; i <= n3h; i++) {
    free_dmatrix(v[i], nrl, nrh, ncl, nch);
  }
  free_new_vvector(v, n3l, n3h, NRUTIL_DPTR2);
}
double ****dmatrix4(unsigned long long n4l, unsigned long long n4h, unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  double ****v = (double ****) new_vvector(n4l, n4h, NRUTIL_DPTR3);
  for(unsigned long long i = n4l; i <= n4h; i++) {
    v[i] = dmatrix3(n3l, n3h, nrl, nrh, ncl, nch);
  }
  return v;
}
void free_dmatrix4(double ****v, unsigned long long n4l, unsigned long long n4h, unsigned long long n3l, unsigned long long n3h, unsigned long long nrl, unsigned long long nrh, unsigned long long ncl, unsigned long long nch) {
  for(unsigned long long i = n4l; i <= n4h; i++) {
    free_dmatrix3(v[i], n3l, n3h, nrl, nrh, ncl, nch);
  }
  free_new_vvector(v, n4l, n4h, NRUTIL_DPTR3);
}
#ifdef _OPENMP
omp_lock_t *ompvector(unsigned long long nl, unsigned long long nh) {
  return ((omp_lock_t *) gvector(nl, nh, sizeof(omp_lock_t)) -nl+NR_END);
}
void free_ompvector(omp_lock_t *v, unsigned long long nl, unsigned long long nh) {
  free_gvector(v+nl-NR_END, nl, nh, sizeof(omp_lock_t));
}
#endif
void *new_vvector(unsigned long long nl, unsigned long long nh, enum alloc_type type) {
  void *v;
  v = NULL;  
  switch(type) {
  case NRUTIL_DPTR:
    v = (double **) gvector(nl, nh, sizeof(double*)) -nl+NR_END;
    break;
  case NRUTIL_UPTR:
    v = (unsigned int **) gvector(nl, nh, sizeof(unsigned int*)) -nl+NR_END;
    break;
  case NRUTIL_DPTR2:
    v = (double ***) gvector(nl, nh, sizeof(double**)) -nl+NR_END;
    break;
  case NRUTIL_NPTR:
    v = (Node **) gvector(nl, nh, sizeof(Node*)) -nl+NR_END;
    break;
  case NRUTIL_NPTR2:
    v = (Node ***) gvector(nl, nh, sizeof(Node**)) -nl+NR_END;
    break;
  case NRUTIL_CPTR:
    v = (char **) gvector(nl, nh, sizeof(char*)) -nl+NR_END;
    break;
  case NRUTIL_DPTR4:
    v = (double *****) gvector(nl, nh, sizeof(double****)) -nl+NR_END;
    break;
  case NRUTIL_TPTR:
    v = (Terminal **) gvector(nl, nh, sizeof(Terminal*)) -nl+NR_END;
    break;
  case NRUTIL_TPTR2:
    v = (Terminal ***) gvector(nl, nh, sizeof(Terminal**)) -nl+NR_END;
    break;
  case NRUTIL_IPTR:
    v = (int **) gvector(nl, nh, sizeof(int*)) -nl+NR_END;
    break;
  case NRUTIL_IPTR2:
    v = (int ***) gvector(nl, nh, sizeof(int**)) -nl+NR_END;
    break;
  case NRUTIL_NPTR3:
    v = (Node ****) gvector(nl, nh, sizeof(Node***)) -nl+NR_END;
    break;
  case NRUTIL_FPTR:
    v = (Factor **) gvector(nl, nh, sizeof(Factor*)) -nl+NR_END;
    break;
  case NRUTIL_FPTR2:
    v = (Factor ***) gvector(nl, nh, sizeof(Factor**)) -nl+NR_END;
    break;
  case NRUTIL_DPTR3:
    v = (double ****) gvector(nl, nh, sizeof(double***)) -nl+NR_END;
    break;
  case NRUTIL_UPTR3:
    v = (unsigned int ****) gvector(nl, nh, sizeof(unsigned int***)) -nl+NR_END;
    break;
  case NRUTIL_UPTR4:
    v = (unsigned int *****) gvector(nl, nh, sizeof(unsigned int****)) -nl+NR_END;
    break;
  case NRUTIL_UPTR2:
    v = (unsigned int ***) gvector(nl, nh, sizeof(unsigned int**)) -nl+NR_END;
    break;
  case NRUTIL_XPTR:
    v = (SNPAuxiliaryInfo **) gvector(nl, nh, sizeof(SNPAuxiliaryInfo*)) -nl+NR_END;
    break;
  case NRUTIL_QPTR:
    v = (QuantileObj **) gvector(nl, nh, sizeof(QuantileObj*)) -nl+NR_END;
    break;
  case NRUTIL_QPTR2:
    v = (QuantileObj ***) gvector(nl, nh, sizeof(QuantileObj**)) -nl+NR_END;
    break;
  case NRUTIL_SPTR:
    v = (LookUpInfo **) gvector(nl, nh, sizeof(LookUpInfo*)) -nl+NR_END;
    break;
  case NRUTIL_SPTR2:
    v = (LookUpInfo ***) gvector(nl, nh, sizeof(LookUpInfo**)) -nl+NR_END;
    break;
  case NRUTIL_VPTR:
    v = (void **) gvector(nl, nh, sizeof(void*)) -nl+NR_END;
    break;
  case NRUTIL_LPTR:
    v = (unsigned long **) gvector(nl, nh, sizeof(unsigned long*)) -nl+NR_END;
    break;
#ifdef _OPENMP
  case NRUTIL_OMPLPTR:
    v = (omp_lock_t **) gvector(nl, nh, sizeof(omp_lock_t*)) -nl+NR_END;
    break;
  case NRUTIL_OMPLPTR2:
    v = (omp_lock_t ***) gvector(nl, nh, sizeof(omp_lock_t**)) -nl+NR_END;
    break;
#endif
  case NRUTIL_LEAFPTR:
    v = (LeafLinkedObj **) gvector(nl, nh, sizeof(LeafLinkedObj*)) -nl+NR_END;
    break;
  case NRUTIL_LEAFPTR2:
    v = (LeafLinkedObj ***) gvector(nl, nh, sizeof(LeafLinkedObj**)) -nl+NR_END;
    break;
  case NRUTIL_SRTLNKPTR:
    v = (SortedLinkedObj **) gvector(nl, nh, sizeof(SortedLinkedObj*)) -nl+NR_END;
    break;
  default:
    v = NULL;
    nrerror("\n  Illegal case in new_vvector().");
    break;
  }
  return v;
}
void free_new_vvector(void *v, unsigned long long nl, unsigned long long nh, enum alloc_type type) {
  switch(type) {
  case NRUTIL_DPTR:
    free_gvector((double**) v +nl-NR_END, nl, nh, sizeof(double*));
    break;
  case NRUTIL_UPTR:
    free_gvector((unsigned int**) v +nl-NR_END, nl, nh, sizeof(unsigned int*));
    break;
  case NRUTIL_DPTR2:
    free_gvector((double***) v +nl-NR_END, nl, nh, sizeof(double**));
    break;
  case NRUTIL_NPTR:
    free_gvector((Node**) v +nl-NR_END, nl, nh, sizeof(Node*));
    break;
  case NRUTIL_NPTR2:
    free_gvector((Node***) v +nl-NR_END, nl, nh, sizeof(Node**));
    break;
  case NRUTIL_CPTR:
    free_gvector((char**) v +nl-NR_END, nl, nh, sizeof(char*));
    break;
  case NRUTIL_DPTR4:
    free_gvector((double*****) v +nl-NR_END, nl, nh, sizeof(double****));
    break;
  case NRUTIL_TPTR:
    free_gvector((Terminal**) v +nl-NR_END, nl, nh, sizeof(Terminal*));
    break;
  case NRUTIL_TPTR2:
    free_gvector((Terminal***) v +nl-NR_END, nl, nh, sizeof(Terminal**));
    break;
  case NRUTIL_IPTR:
    free_gvector((int**) v +nl-NR_END, nl, nh, sizeof(int*));
    break;
  case NRUTIL_IPTR2:
    free_gvector((int***) v +nl-NR_END, nl, nh, sizeof(int**));
    break;
  case NRUTIL_NPTR3:
    free_gvector((Node****) v +nl-NR_END, nl, nh, sizeof(Node***));
    break;
  case NRUTIL_FPTR:
    free_gvector((Factor**) v +nl-NR_END, nl, nh, sizeof(Factor*));
    break;
  case NRUTIL_FPTR2:
    free_gvector((Factor***) v +nl-NR_END, nl, nh, sizeof(Factor**));
    break;
  case NRUTIL_DPTR3:
    free_gvector((double****) v +nl-NR_END, nl, nh, sizeof(double***));
    break;
  case NRUTIL_UPTR3:
    free_gvector((unsigned int****) v +nl-NR_END, nl, nh, sizeof(unsigned int***));
    break;
  case NRUTIL_UPTR4:
    free_gvector((unsigned int*****) v +nl-NR_END, nl, nh, sizeof(unsigned int****));
    break;
  case NRUTIL_UPTR2:
    free_gvector((unsigned int***) v +nl-NR_END, nl, nh, sizeof(unsigned int**));
    break;
  case NRUTIL_XPTR:
    free_gvector((SNPAuxiliaryInfo **) v +nl-NR_END, nl, nh, sizeof(SNPAuxiliaryInfo*));
    break;
  case NRUTIL_QPTR:
    free_gvector((QuantileObj **) v +nl-NR_END, nl, nh, sizeof(QuantileObj*));
    break;
  case NRUTIL_QPTR2:
    free_gvector((QuantileObj ***) v +nl-NR_END, nl, nh, sizeof(QuantileObj**));
    break;
  case NRUTIL_SPTR:
    free_gvector((LookUpInfo **) v +nl-NR_END, nl, nh, sizeof(LookUpInfo*));
    break;
  case NRUTIL_SPTR2:
    free_gvector((LookUpInfo ***) v +nl-NR_END, nl, nh, sizeof(LookUpInfo**));
    break;
  case NRUTIL_VPTR:
    free_gvector((void**) v +nl-NR_END, nl, nh, sizeof(void*));
    break;
  case NRUTIL_LPTR:
    free_gvector((unsigned long**) v +nl-NR_END, nl, nh, sizeof(unsigned long*));
    break;
#ifdef _OPENMP
  case NRUTIL_OMPLPTR:
    free_gvector((omp_lock_t**) v +nl-NR_END, nl, nh, sizeof(omp_lock_t*));
    break;
  case NRUTIL_OMPLPTR2:
    free_gvector((omp_lock_t***) v +nl-NR_END, nl, nh, sizeof(omp_lock_t**));
    break;
#endif
  case NRUTIL_LEAFPTR:
    free_gvector((LeafLinkedObj**) v +nl-NR_END, nl, nh, sizeof(LeafLinkedObj*));
    break;
  case NRUTIL_LEAFPTR2:
    free_gvector((LeafLinkedObj***) v +nl-NR_END, nl, nh, sizeof(LeafLinkedObj**));
    break;
  case NRUTIL_SRTLNKPTR:
    free_gvector((SortedLinkedObj**) v +nl-NR_END, nl, nh, sizeof(SortedLinkedObj*));
    break;
  default:
    nrerror("\n  Illegal case in free_new_vvector().");
    break;
  }
}
#undef FREE_ARG
void nrCopyMatrix(unsigned int **new, unsigned int **old, unsigned int nrow, unsigned int ncol) {
  unsigned int i,j;
  for (i = 1; i <= nrow; i++) {
    for (j = 1; j <= ncol; j++) {
      new[i][j] = old[i][j];
    }
  }
}
void nrCopyVector(char *new, char *old, unsigned int ncol) {
  unsigned int j;
  for (j = 1; j <= ncol; j++) {
    new[j] = old[j];
  }
}
void testEndianness(void) {
  unsigned int     test = 0x12345678;
  unsigned int *testPtr = & test;
  RF_nativePrint("\nTest of Endianness:  ");
  RF_nativePrint("%2x %2x %2x %2x \n",
           *((char *) testPtr),
           *((char *) testPtr + 1),
           *((char *) testPtr + 2),
           *((char *) testPtr + 3));
}
