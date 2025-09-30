
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "svdUtil.h"
#include "nrutil.h"
#include "error.h"
void svdcmp(double **aorg, int m, int n, double ***uptr, double **wptr, double ***vptr) {
  double **a, *w, **v;
  double pythag(double a, double b);
  int flag, i, its, j, jj, k, l, nm;
  double anorm, c, f, g, h, s, scale, x, y, z, *rv1;
  double sgf;
  *uptr = a = matrixCopy(aorg, m, n);
  *wptr = w = dvector(1, n);
  *vptr = v = dmatrix(1, n, 1, n);
  rv1 = dvector(1, n);
  g = scale = anorm = 0.0; 
  for (i = 1; i <= n; i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
    if (i <= m) {
      for (k = i; k <= m; k++) scale += fabs(a[k][i]);
      if (scale) {
        for (k = i; k <= m; k++) {
          a[k][i] /= scale;
          s += a[k][i] * a[k][i];
        }
        f = a[i][i];
        if (f >= 0.0) {
          g = -sqrt(s);
        }
        else {
          g = sqrt(s);
        }
        h = (f * g) - s;
        a[i][i] = f - g;
        for (j = l; j <= n; j++) {
          for (s = 0.0, k = i; k <= m; k++) s += a[k][i] * a[k][j];
          f = s / h;
          for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
        }
        for (k = i; k <= m; k++) a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if ((i <= m) && (i != n)) {
      for (k = l; k <= n; k++) scale += fabs(a[i][k]);
      if (scale) {
        for (k = l; k <= n; k++) {
          a[i][k] /= scale;
          s += a[i][k] * a[i][k];
        }
        f = a[i][l];
        if (f >= 0.0) {
          g = -sqrt(s);
        }
        else {
          g = sqrt(s);
        }
        h = (f * g) - s;
        a[i][l] = f - g;
        for (k = l; k <= n; k++) rv1[k] = a[i][k]/h;
        for (j = l; j <= m; j++) {
          for (s = 0.0,k = l; k <= n; k++) s += a[j][k]*a[i][k];
          for (k = l; k <= n; k++) a[j][k] += s*rv1[k];
        }
        for (k = l; k <= n; k++) a[i][k] *= scale;
      }
    }
    anorm = (anorm > (fabs(w[i]) + fabs(rv1[i]))) ? anorm : (fabs(w[i]) + fabs(rv1[i]));
  }
  for (i = n; i >= 1; i--) { 
      if (i < n) {
        if (g) {
          for (j = l; j <= n; j++) 
            v[j][i] = (a[i][j]/a[i][l]) / g;
          for (j = l; j <= n; j++) {
            for (s = 0.0, k = l; k <= n; k++) s += a[i][k]*v[k][j];
            for (k = l; k <= n; k++) v[k][j] += s*v[k][i];
          }
        }
        for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
      }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }
  for (i = (m < n) ? m : n; i >= 1; i--) { 
    l = i+1;
    g = w[i];
    for (j = l; j <= n; j++) a[i][j] = 0.0;
    if (g) {
      g = 1.0/g;
      for (j = l; j <= n; j++) {
        for (s = 0.0,k = l; k <= m; k++) s += a[k][i]*a[k][j];
        f = (s/a[i][i])*g;
        for (k = i; k <= m; k++) a[k][j] += f*a[k][i];
      }
      for (j = i; j <= m; j++) a[j][i] *= g;
    } else for (j = i; j <= m; j++) a[j][i] = 0.0;
    ++a[i][i];
  }
  for (k = n; k >= 1; k--) {
    for (its = 1; its <= 30; its++) {
      flag = 1;
      for (l = k; l >= 1; l--) { 
          nm = l - 1;
          if ((double) (fabs(rv1[l]) + anorm) == anorm) {
          flag = 0;
          break;
        }
        if ((double) (fabs(w[nm]) + anorm) == anorm) break;
      }
      if (flag) {
        c = 0.0; 
        s = 1.0;
        for (i = l; i <= k; i++) {
          f = s*rv1[i];
          rv1[i] = c*rv1[i];
          if ((double) (fabs(f) + anorm) == anorm) break;
          g = w[i];
          h = pythag(f,g);
          w[i] = h;
          h = 1.0/h;
          c = g*h;
          s = -f*h;
          for (j = 1; j <= m; j++) {
            y = a[j][nm];
            z = a[j][i];
            a[j][nm] = y*c+z*s;
            a[j][i] = z*c-y*s;
          }
        }
      }
      z = w[k];
      if (l == k) {
        if (z < 0.0) { 
          w[k] = -z;
          for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == 30) nrerror("no convergence in 30 SVD iterations");
      x = w[l]; 
      nm = k-1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ( ((y-z) * (y+z)) + ((g-h) * (g+h)) ) / (2.0 * h * y);
      g = pythag(f,1.0);
      if (f >= 0.0) {
        sgf = fabs(g);
      }
      else {
        sgf = -fabs(g);
      }
      f = ( ((x-z) * (x+z)) + h * ( (y / ( f + sgf)) - h) ) / x;      
      c = s = 1.0; 
      for (j = l; j <= nm; j++) {
        i = j+1;
        g = rv1[i];
        y = w[i];
        h = s*g;
        g = c*g;
        z = pythag(f,h);
        rv1[j] = z;
        c = f/z;
        s = h/z;
        f = x*c+g*s;
        g  =  g*c-x*s;
        h = y*s;
        y *= c;
        for (jj = 1; jj <= n; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }
        z = pythag(f,h);
        w[j] = z; 
        if (z) {
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }
        f = c * g + s * y;
        x = c * y - s * g;
        for (jj = 1; jj <= m; jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y * c + z * s;
          a[jj][i] = z * c - y * s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  free_dvector(rv1, 1, n);
}
char svdchk(double **a, uint m, uint n, double **u, double *w, double **v) {
  double **atest;
  uint k, j, i;
  double **tmp;
  char result;
  atest = dmatrix (1, m, 1, n);
  tmp = dmatrix(1, m, 1, n);
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++) {
      tmp[i][j] = u[i][j] * w[j];
    }
  }
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++) {
      atest[i][j] = 0.0;
      for (k = 1; k <= n; k++) {
        atest[i][j] += tmp[i][k] * v[j][k];
      }
    }
  }
  free_dmatrix(tmp, 1, m, 1, n);
  RF_nativePrint("\n");
  RF_nativePrint("\n Original [A] of dim m x n :");
  matrixPrint(a, m, n);
  RF_nativePrint("\n");
  RF_nativePrint("\n Recovered [A] of dim m x n :");
  matrixPrint(atest, m, n);
  result = TRUE;
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++) {
      if (fabs(atest[i][j] - a[i][j]) > EPSILON) {
        result = FALSE;
       }
    }
  }
  RF_nativePrint("\n");
  if (result) {
    RF_nativePrint("\n Original [A] == Recovered [A] ? : TRUE");
  }
  else {
    RF_nativePrint("\n Original [A] == Recovered [A] ? : FALSE");
  }
  free_dmatrix(atest, 1, m, 1, n);
  return result;
}
double **svdinv(double **u, double *w, double **v, uint m, uint n, uint singularity) {
  double **aplus, **wplus, **utrans, **tmp;
  uint i, j, k;
  wplus = dmatrix(1, n, 1, n);
  k = 0;
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      if (i != j) {
        wplus[i][j] = 0.0;
      }
      else {
        if (fabs(w[i]) > SVD_EPS) {
          wplus[i][j] = 1.0 / w[i];
          k++;
        }
        else {
          wplus[i][j] = 0.0;
        }
      }
    }
  }
  if ((k >= singularity) && k > 1) {
    tmp = matrixMult(v, wplus, n, n, n);
    utrans = matrixTrans(u, m, n);
    aplus = matrixMult(tmp, utrans, n, n, m);
    free_dmatrix(tmp, 1, n, 1, n);
    free_dmatrix(utrans, 1, n, 1, m);
  }
  else {
    aplus = NULL;
  }
  free_dmatrix(wplus, 1, n, 1, n);
  return aplus;
}
void free_svdcmp(double **a, int m, int n, double **u, double *w, double **v) {
  if (a != NULL) {
    free_dmatrix(a, 1, m, 1, n);
    a = NULL;
  }
  free_dmatrix(u, 1, m, 1, n);
  u = NULL;
  free_dvector(w, 1, n);
  w = NULL;
  free_dmatrix(v, 1, n, 1, n);
  v = NULL;
}
void svbksb(double **u, double *w, double **v, uint m, uint n, double *b, double *x) {
  uint jj, j, i;
  double s, *tmp;
  tmp = dvector(1, n);
  for (j=1; j <= n; j++) {
    s = 0.0;
    if (w[j]) {
      for (i=1; i <= m; i++) s += u[i][j] * b[i];
      s = s / w[j]; 
    }
    tmp[j]=s;
  }
  for (j = 1; j <= n; j++) {
    s = 0.0;
    for (jj = 1; jj <= n; jj++) s += v[j][jj] * tmp[jj];
    x[j] = s;
  }
  free_dvector(tmp, 1, n);
}
double **matrixCopy(double **a, uint m, uint n) {
  double **acopy;
  uint i, j;
  acopy = dmatrix(1, m, 1, n);
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++) {
      acopy[i][j] = a[i][j];
    }
  }
  return acopy;
}
double **matrixTrans(double **a, uint m, uint n) {
  uint i, j;
  double **atrans;
  atrans = dmatrix(1, n, 1, m);
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++) {
      atrans[j][i] = a[i][j];
    }
  }
  return atrans;
}
double **matrixMult(double **a, double **b, uint m, uint n, uint p) {
  double **c;
  uint i, j, k;
  c = dmatrix(1, m, 1, p);
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= p; j++) {
      c[i][j] = 0.0;
      for (k = 1; k <= n; k++) {
        c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return c;
}
void matrixPrint(double **x, uint m, uint n) {
  uint i, j;
  for (i = 1; i <= m; i++) {
    RF_nativePrint("\n");
    for (j = 1; j <= n; j++) {
      RF_nativePrint("  %10.8e", x[i][j]);
    }
  }
}
double pythag(double a, double b) {
  double absa, absb;
  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb)
    return absa * sqrt(1.0 + ((absb/absa) * (absb/absa)));
  else
    return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + ((absa/absb) * (absa/absb))));
}
void harness(void) {
  uint M, N;
  M = 3;
  N = 4;
  double **a, **aplus, **gident, **gident2, **uident, **uident2;
  double  *w;
  double **u, **utrans;
  double **v, **vtrans, **vident, **vident2;
  a    = dmatrix(1, M, 1, N);
  a[1][1] =  3;
  a[1][2] =  0;
  a[1][3] =  1;
  a[1][4] =  3;
  a[2][1] =  4;
  a[2][2] =  5;
  a[2][3] =  6;
  a[2][4] =  7;
  a[3][1] =  0;
  a[3][2] =  1;
  a[3][3] =  2;
  a[3][4] = -2;
  svdcmp(a, M, N, &u, &w, &v);
  RF_nativePrint("\n");  
  RF_nativePrint("\n Original [A] of dim m x n :");
  matrixPrint(a, M, N);
  svdchk(a, M, N, u, w, v);
  aplus = svdinv(u, w, v, M, N, 2);
  gident = matrixMult(a, aplus, M, N, M);
  RF_nativePrint("\n");
  RF_nativePrint("\n Check of [A] x [A+] of dim m x m :");
  matrixPrint(gident, M, M);
  gident2 = matrixMult(aplus, a, N, M, N);
  RF_nativePrint("\n");
  RF_nativePrint("\n Check of [A+] x [A] of dim n x n :");
  matrixPrint(gident2, N, N);
  utrans = matrixTrans(u, M, N);
  uident = matrixMult(u, utrans, M, N, M);
  RF_nativePrint("\n");
  RF_nativePrint("\n Check of [U] x [U]^t of dim m x m :");
  matrixPrint(uident, M, M);
  uident2 = matrixMult(utrans, u, N, M, N);
  RF_nativePrint("\n");
  RF_nativePrint("\n Check of [U]^t x [U] of dim n x n :");
  matrixPrint(uident2, N, N);
  vtrans = matrixTrans(v, N, N);
  vident = matrixMult(v, vtrans, N, N, N);
  RF_nativePrint("\n");
  RF_nativePrint("\n Check of [V] x [V]^t of dim n x n :");
  vident2 = matrixMult(vtrans, v, N, N, N);
  RF_nativePrint("\n");
  RF_nativePrint("\n Check of [V]^t x [V] of dim n x n :");
  matrixPrint(vident2, N, N);
  free_dmatrix(vident2, 1, N, 1, N);
  free_dmatrix(vident, 1, N, 1, N);
  free_dmatrix(vtrans, 1, N, 1, N);
  free_dmatrix(uident2, 1, N, 1, N);
  free_dmatrix(uident, 1, M, 1, M);
  free_dmatrix(utrans, 1, N, 1, M);
  free_dmatrix(gident2, 1, N, 1, N);
  free_dmatrix(gident, 1, M, 1, M);  
  free_dmatrix(aplus, 1, N, 1, M);  
  free_svdcmp(a, M, N, u, w, v);
}
