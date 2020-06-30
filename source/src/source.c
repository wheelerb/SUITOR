
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


#define binary char
#define NMF_STOPCONV 40
#define NMF_CHECKINTERVAL 10
#define NMF_MAXITER 2000
#define CHECK_MEM(obj) if (obj == NULL) {Rprintf("ERROR: allocating memory \n"); error("1");}
#define NUMERICZERO 1e-200
#define LARGEDOUBLE 1.0e200
#define MINUSINFINITY -1.0e200
#define DOUBLE_MISS -9999.0e100
#define MAX( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )


#define IARG_NR 0
#define IARG_NC 1
#define IARG_RANK 2
#define IARG_KFOLD 3
#define IARG_EM_MAXITER 4
#define IARG_MAX_RANK 5
#define IARG_PRINT 6
#define IARG_NUM_K 7

#define DARG_UPDATE_EPS 0
#define DARG_EM_EPS 1
#define DARG_EM_MINVALUE 2


void C_call_nmf(double *mat, int *iargs, double *dargs, double *retW, double *retH);
void C_call_suitor(double *mat, int *iargs, double *dargs, int *rankVec, int *kVec, 
                   double *ret_train, double *ret_test);

static const R_CMethodDef callMethods[] = {
  {"C_call_nmf", (DL_FUNC)&C_call_nmf, 5},
  {"C_call_suitor", (DL_FUNC)&C_call_suitor, 7},
  {NULL, NULL, 0}
};


struct nmf_struct {
  
  double *x;       /* input mat as a vector (columns stacked) */
  double *logx;    /* log(x) if x > NUMERICZERO, 0 otherwise */
  double *x_k;
  double *x_k_hat;
  double *x_k0;
  double *h, *h0;
  double *w, *w0;
  int xnr, xnc, hnr, hnc, wnr, wnc; /* Number of rows cols for x, h, w */
  int xlen, wlen, hlen;  /* xnr*xnc, etc */

  int rank;      /* desired rank */
  int kfold;     /* number of folds */
  int *rankVec;
  int *kVec;
  int nk;           /* length of kVec and rankVec */
  double *runifVec; /* Vector of uniform(0,1) value for W and H */
  int maxRank;      /* maximum value in rankVec */
  
  binary *cons;      /* Connectivity matrix */
  int cons_len;      /* hnc*(hnc-1)/2 */
  int *h_max_index;  /* Index for max value in each column of h */
  int cons_inc;      /* counting identical cons */
  double update_eps; /* Min value for w, h mats */
  double EM_eps;     /* Stop tol in EM */
  int EM_maxiter;
  binary *idxMat;  
  int *idxInt;        /* Vector of integer positions in idxMat */  
  double EM_minvalue;  /* replace small values with min_value in input.k in EM alg */
  int EM_iter;         /* Number of iterations for convergence */
  double delta_denom;       
  double train_adj;
  double loglike;
  int print;
  double *train_errors;
  double *test_errors;

  double *tmpVecN, *tmpVecN2; /* scratch vectors of length xlen */
  double *tmpVecM;            /* scratch vector of length max(xnr, xnc) */

};
typedef struct nmf_struct NMFSTR;
 

static void print_dVec(vec, n, name, by)
double *vec;
int n, by;
char name[10];
{
  int i, j=0;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %g ", vec[i]);
    j++;
    if (by && (j >= by)) {
      Rprintf("\n");
      j = 0;
    }
  }
  Rprintf("\n");
}

static void print_dMat(vec, nr, nc, name)
double *vec;
int nr, nc;
char name[10];
{
  int i, j;
  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      Rprintf(" %g ", vec[i + j*nr]);
    }
    Rprintf("\n");
  }
  
}

static void print_bMat(vec, nr, nc, name)
binary *vec;
int nr, nc;
char name[10];
{
  int i, j, k;
  binary tmp;

  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      tmp = vec[i + j*nr];
      if (!tmp) {
        k = 0;
      } else {
        k = 1;
      } 
      Rprintf(" %d ", k);
    }
    Rprintf("\n");
  }
  
}


static void print_iVec(vec, n, name)
int *vec;
int n;
char name[10];
{
  int i;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %d ", vec[i]);
  }
  Rprintf("\n");
}



static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  if (n < 1) error("n < 1 in dVec_alloc");
  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} 

static int * iVec_alloc(n, initFlag, initVal)
int n, initFlag, initVal;
{
  int i, *ret, *p;

  if (n < 1) error("n < 1 in iVec_alloc");
  ret = (int *) malloc(n*sizeof(int));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} 


static char * cVec_alloc(n, initFlag, initVal)
int n, initFlag;
char initVal;
{
  char *ret, *pret;
  int i;

  if (n < 1) error("n < 1 in cVec_alloc");
  ret = (char *) malloc(n*sizeof(char));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, pret=ret; i<n; i++, pret++) *pret = initVal;
  }
  
  return(ret);

} 

static void copy_dVec(v, n, ret)
double *v, *ret;
int n;
{
  int i;
  double *p1, *p2;

  for (i=0, p1=v, p2=ret; i<n; i++, p1++, p2++) *p2 = *p1;

}

static void matrixMult(m1, m1_nr, m1_nc, m2, m2_nc, ret)
double *m1, *m2, *ret;
int m1_nr, m1_nc, m2_nc;
{
  /*  all matrices are stacked columns, same for ret */

  int i, j, k, ii, colind;
  double sum;

  for (i=0; i<m1_nr; i++) {
    for (j=0; j<m2_nc; j++) {
      sum    = 0.0;
      colind = j*m1_nc;  
      for (k=0; k<m1_nc; k++) {
        sum += m1[i + m1_nr*k]*m2[colind];
        colind++;
      }
      ret[i + j*m1_nr] = sum;
    }
  }

}

static int which_max(v, n) 
double *v;
int n;
{

  int i, maxi;
  double maxd, tmp;

  maxi = 0;
  maxd = v[0];

  if (n > 1) {
    for (i=1; i<n; i++) {
      tmp = v[i];
      if (tmp > maxd) {
        maxd = tmp;
        maxi = i;
      }
    }
  }

  return(maxi);
}

static void which_max_cols(mat, nr, nc, ret) 
double *mat; /* stacked columns */
int nr, nc, *ret;
{
  int i;

  for (i=0; i<nc; i++) {
    ret[i] = which_max(&mat[i*nr], nr); 
  }

}

static int update_cons(cons, cons_len, index, nc)
binary *cons;  /* cons is updated in this function */
int cons_len, *index, nc; 
{
  /*cons  <- outer(index, index, function(x,y) ifelse(x==y, 1,0))*/
  int i, j, indi, indj, cons_index, anyChange=0;
  binary isEq;

  /* Loop over upper triangle of outer (not including diagonal) */
  cons_index = 0;
  for (i=0; i<nc-1; i++) {
    indi = index[i];
    for (j=i+1; j < nc; j++) {
      if (indi == index[j]) {
        isEq = 1;
      } else {
        isEq = 0;
      }
      if (cons[cons_index] != isEq) anyChange = 1;
      if (cons_index >= cons_len) error("INTERNAL CODING ERROR in update_cons");
      cons[cons_index] = isEq;
      cons_index++;
      
    }
  }

  return(anyChange);

}

/* cons and inc get updated */
static int myStopFun(cons, cons_len, inc, hMat, hnr, hnc, maxIndex) 
binary *cons;
int cons_len, *inc, hnr, hnc, *maxIndex;
double *hMat;
{
  int anyChange, ret=0;

  /* Get the maximum in each column of H */
  which_max_cols(hMat, hnr, hnc, maxIndex);

  /* Update connectivity matrix */
  anyChange = update_cons(cons, cons_len, maxIndex, hnc);

  if (anyChange) {
    *inc = 0;
  } else {
    *inc = *inc + 1;
  }

  if (*inc > NMF_STOPCONV) ret = 1;

  return(ret);

}

static void get_runif_vec(v, n, a, b)
double *v, a, b;
int n;
{
  int i;

  for (i=0; i<n; i++) v[i] = runif(a, b);

}

/* Return sum of idxMat (used later) */
static int get_idxMat0(idx, nr, nc, k, kfold) 
binary *idx;
int nr, nc, k, kfold;
{
  int i, val, j, offset, ret=0;

  for (i=0; i<nr*nc; i++) idx[i] = 0;
  for (i=1; i<=nc; i++) {
    offset = (i-1)*nr;
    val = (i + k - 1) % kfold;
    if (!val) val += kfold;

    while (1) {
      if (val <= nr) {   
        idx[val-1+offset] = 1;
        ret++;
      } else {
        break;
      }
      val += kfold;
    } 
  }

  return(ret);
}

static double loglike(input, log_input, input_k_hat, idxMat, delta_denom, len, which_idx)
double *input, *log_input, *input_k_hat, delta_denom;
int len;
binary *idxMat, which_idx;
{
   int i;
   double ret=0.0, inputi, lik, inputkhati, tmp, logDD;
  
   logDD = log(delta_denom);
   for (i=0; i<len; i++) {  
     if (idxMat[i] == which_idx) {
       inputi     = input[i];
       inputkhati = input_k_hat[i];
       if (inputi > NUMERICZERO) {       
         if (inputkhati < delta_denom) {
           tmp = logDD;
         } else {
           tmp = log(inputkhati);
         }
         lik  = inputi*(log_input[i] - tmp);
         ret += lik - inputi + inputkhati;
       } else {
         ret += inputkhati;
       }     
     }
   }
   
   return(ret);

}

static void replaceZeroVal_dVec(v, n, minval)
double *v, minval;
int n;
{
  int i;

  for (i=0; i<n; i++) {
    if (v[i] < minval) v[i] = minval;
  }

}

static double get_delta_denom(v, n) 
double *v;
int n;
{
  int i;
  double minv=LARGEDOUBLE, tmp;

  for (i=0; i<n; i++) {
    tmp = v[i];
    if ((tmp > NUMERICZERO) && (tmp < minv)) minv = tmp;
  }
  minv *= 0.5;

  return(minv);
}

static double get_max(v, n)
double *v;
int n;
{
  int i;
  double maxv=MINUSINFINITY, tmp;
  
  for (i=0; i<n; i++) {
    tmp = v[i];
    if (tmp > maxv) maxv = tmp; 
  }

  return(maxv);

}

static void my_update_H(pV, w, h, wnr, wnc, hnc)
double *pV, *w, *h;
int wnr, wnc, hnc;
{
   int n, r, p, ncterms=0, nbterms=0, vr, iH, jH, u, k, ii;
   double *res, *pW, *pH, *p_res, *sumW, *pWH, tmp_res, w_sum, wh_term;

   /* retrieve dimensions from W and H */
   n   = wnr;
   r   = wnc;
   p   = hnc;
   res = h;

   /* get number of non-fixed terms */
   vr = r - ncterms;

   /* define internal pointers */
   pW    = w;
   pH    = h;
   p_res = res;

   /* allocate internal memory */
   sumW = dVec_alloc(r, 0, 0.0); /* will store column sums of W */
   pWH  = dVec_alloc(n, 0, 0.0); /* will store the currently used column of WH */

   /* Get column sums of w */
   ii = 0;
   for (u=0; u<r; u++) {
     w_sum = 0.0;
     for (k=0; k<n; k++) {
       w_sum += w[ii];
       ii++;
     }
     sumW[u] = w_sum;
   }


       /* Compute update of H column by column */
	for (jH=0; jH < p; ++jH){

		for (iH=0; iH < vr; ++iH){ /* compute value for H_ij (non-fixed terms only) */

			/* initialize values */
			tmp_res = 0.0;
			w_sum   = sumW[iH];
			/*if( jH == 0 ) w_sum = 0.0;*/

			/* compute cross-product w_.i by (v/wh)_.j */
			for(u=0; u<n; u++){

				/* The jH-th column of WH is used to compute all elements of H_.j
				    => compute once and store the result for using for the next rows */
				wh_term = pWH[u];
				if( iH == 0 ){
					wh_term = 0.0;
					for (k=0; k<r; k++){
						wh_term += pW[u + k*n] * pH[k + jH*r];
					}
					wh_term = pV[u + jH*n] / wh_term;
					pWH[u] = wh_term;
				}

				tmp_res +=  pW[u + iH*n] * wh_term;

				/* compute sum of iH-th column of W (done only once)
				if( jH == 0 ) w_sum += pW[u + iH*n];
                            */
			}

			/* multiplicative update */
			p_res[iH + jH*r] = pH[iH + jH*r] * tmp_res / w_sum;
		}
	}

  free(sumW);
  free(pWH);

}

static void my_update_W(pV, w, h, wnr, wnc, hnc)
double *pV, *w, *h;
int wnr, wnc, hnc;
{

	double *res, *pW, *pH, *pWH, *p_res, *sumH, wh_term, h_sum, tmp_res;
       int n, r, p, iW, u, jW, k;

	/* retrieve dimensions */
	n = wnr;
	r = wnc;
	p = hnc;

       res = w;

	/* define internal pointers */
	pW    = w;
	pH    = h;
	p_res = res;

	/* allocate internal memory */
	sumH = dVec_alloc(r, 0, 0.0); /* will store the row sums of H */
	pWH =  dVec_alloc(p, 0, 0.0); /* will store currently used row of WH */

       /* row sums of h */
       for (jW=0; jW < r; jW++) {
          h_sum = 0.0;
          for(u=0; u<p; u++) h_sum += pH[jW + u*r];     
          sumH[jW] = h_sum;
       }
      

	/* Compute update of W row by row */
	for(iW=0; iW < n; iW++){

		for (jW=0; jW < r; jW++){ /* compute value for W_ij */

			/* initialize values */
			tmp_res = 0.0;
			h_sum   = sumH[jW];
			/*if( iW == 0 ) h_sum = 0.0;*/

			/* compute cross-product (v/wh)_i. by h_j. */
			for(u=0; u<p; u++){

				/* The iW-th row of WH is used to compute all elements of W_i.
				 => compute once and store the result for using for the next columns */
				if( jW == 0 ){
				       wh_term = 0.0;
					for (k=0; k<r; k++){
						wh_term += pW[iW + k*n] * pH[k + u*r];
					}
					wh_term = pV[iW + u*n] / wh_term;
					pWH[u] = wh_term;
				}

				tmp_res +=  pH[jW + u*r] * pWH[u];

				/* compute sum of jW-th row of H (done only once) */
				/*if( iW == 0 ) h_sum += pH[jW + u*r];*/
			}

			/* multiplicative update */
			p_res[iW + jW*n] = pW[iW + jW*n] * tmp_res / h_sum;
		}

	}

	free(sumH);
       free(pWH);
}

static void my_update_brunet(i, x, w, h, str)
int i;
double *x, *w, *h;
NMFSTR *str;
{
  int j, wnr=str->wnr, wnc=str->wnc, hnc=str->hnc;
  double eps=str->update_eps;

  /* standard divergence-reducing NMF update for H */
  my_update_H(x, w, h, wnr, wnc, hnc);

  /* standard divergence-reducing NMF update for W */
  my_update_W(x, w, h, wnr, wnc, hnc);

  if (i % NMF_CHECKINTERVAL == 0) {
    for (j=0; j<str->hlen; j++) {
      if (h[j] < eps) h[j] = eps;
    }
    for (j=0; j<str->wlen; j++) {
      if (w[j] < eps) w[j] = eps;
    }
  }

}

static void get_init_WH(runifVec, w, h, wlen, hlen, maxval)
double *runifVec, *w, *h, maxval;
int wlen, hlen;
{
  int i, j=0;

  /* initialize W first */
  for (i=0; i<wlen; i++) {
    w[i] = maxval*runifVec[j];
    j++;
  }

  for (i=0; i<hlen; i++) {
    h[i] = maxval*runifVec[j];
    j++;
  }

}

static void call_nmf(mat, str)
double *mat;
NMFSTR *str;
{
  /* seed must be set before .C call */

  double maxval, *h, *w;
  int inc, i, conv, len, cons_len, hnr, hnc, *maxIndex;
  binary *cons;

  len      = str->xlen;
  cons     = str->cons;
  cons_len = str->cons_len;
  inc      = 0;
  i        = 0;
  hnr      = str->hnr;
  hnc      = str->hnc;
  maxIndex = str->h_max_index;
  w        = str->w;
  h        = str->h;

  /* get max value of mat */
  maxval = get_max(mat, len);

  get_init_WH(str->runifVec, w, h, str->wlen, str->hlen, maxval);

  while (1) {
    i++;

    /* update w and h */
    my_update_brunet(i, mat, w, h, str);

    /* test convergence only every 10 iterations */
    if (i % NMF_CHECKINTERVAL == 0) {
      conv = myStopFun(cons, cons_len, &inc, h, hnr, hnc, maxIndex);
    }

    if (conv || (i >= NMF_MAXITER)) break;
  }

}

static double get_train_adj(nrnc, idxMat, idxLen)
int nrnc, idxLen;
binary *idxMat;
{
  int i, sum=0;
  double ret;

  for (i=0; i<idxLen; i++) {
    if (idxMat) sum++;
  }
  ret = (double) nrnc - sum;

  return(ret);
}

static void update_x_k(x_k, x_k_hat, idxMat, n)
double *x_k, *x_k_hat;
binary *idxMat;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (idxMat[i]) x_k[i] = x_k_hat[i];
  }
}

/* Sorting function, x is changed, sort from small to large.
   cvec is commented out */
static void quicksort(int start, int stop, double *x)
{
 int i, j, k;
 double temp, median;
 /*int tempd;*/


  while (start < stop) {
    /*
    ** first-- if the list is short, do an ordinary insertion sort
    */
    if ((stop-start)<11) {
	for (i=start+1; i<=stop; i++) {
	    temp = x[i];
	    /*tempd= cvec[i];*/
	    j=i-1;

	    while (j>=start && (x[j]>temp)) {
		x[j+1] = x[j];
		/*cvec[j+1] = cvec[j];*/
		j--;
		}
	    x[j+1] = temp;
	    /*cvec[j+1]  = tempd;*/
	    }
	return;
	}

    /*
    ** list is longer -- split it into two
    **  I use the median of 3 values as the split point
    */
    i=start;
    j=stop;
    k = (start + stop)/2;

    median = x[k];
    if (x[i] >= x[k]) {      /* one of j or k is smallest */
	if (x[j] > x[k]) {   /* k is smallest */
	    if (x[i] > x[j])  median = x[j];
	    else median= x[i];
	    }
	}
    else {
	if (x[j] < x[k]) {
	    if (x[i] > x[j]) median = x[i];
	    else median = x[j];
	    }
	}

    /* 
    **  Now actually do the partitioning 
    **   Because we must have at least one element >= median, "i"
    **   will never run over the end of the array.  Similar logic
    **   applies to j.
    ** A note on the use of "<" rather than "<=".  If a list has lots
    **   of identical elements, e.g. 80/100 are "3.5", then we will
    **   often go to the swap step with x[i]=x[j]=median.  But we will
    **   get the pointers i and j to meet approximately in the middle of
    **   the list, and that is THE important condition for speed in a
    **   quicksort.
    **   
    */
    while (i<j) {
	/*
	** top pointer down till it points at something too large
	*/
	while (x[i] < median) i++;

	/*
	** bottom pointer up until it points at something too small
	*/
	while(x[j] > median) j--;

	if (i<j) {
	    if (x[i] > x[j]) {  /* swap */
		temp = x[i];
		x[i] = x[j];
		x[j] = temp;
		/*tempd= cvec[i];   cvec[i] =cvec[j];  cvec[j] =tempd;*/
		}
	    i++; j--;
	    }
	}

    /*
    ** The while() step helps if there are lots of ties.  It will break
    **  the list into 3 parts: < median, ==median, >=median, of which only
    **  the top and bottom ones need further attention.
    ** The ">=" is needed because i may be  == to j
    */
    while (x[i] >= median && i>start) i--;
    while (x[j] <= median && j<stop ) j++;

    /*
    ** list has been split, now do a recursive call
    **   always recur on the shorter list, as this keeps the total
    **       depth of nested calls to less than log_base2(n).
    */
    if ((i-start) < (stop-j)) { /* top list is shorter */
	if ((i-start)>0) quicksort(start,i, x);
	start =j; 
	}

    else {    /* bottom list is shorter */
	if ((stop -j)>0) quicksort(j,stop, x);
	stop=i; 
	}
     }
}

static double get_median(v, n)
double *v;
int n;
{
  int m;
  double ret;

  /* sort v */
  quicksort(0, n-1, v);

  /* Adjust index since c begins at 0 */

  if (n % 2) {
    /* odd number of values */
    m   = (n + 1)/2 - 1;
    ret = v[m];
  } else {
    m   = n/2 - 1;
    ret = 0.5*(v[m] + v[m+1]);
  }

  return(ret);

}	

static void get_init_x_k(x, x_k, idxMat, nr, nc, tmpvec, minvalue)
double *x, *x_k, *tmpvec, minvalue;
binary *idxMat;
int nr, nc;
{
  int i, j, ii, len, vii;
  double med;

  len = nr*nc;

  for (i=0; i<len; i++) x_k[i] = x[i];

  for (i=0; i<nr; i++) {
    vii = 0;
    for (j=0; j<nc; j++) {
      ii = i + j*nr;
      if (!idxMat[ii]) {
        tmpvec[vii] = x_k[ii];
        vii++;
      }
    }
    /* Get median of tmpvec */
    med = get_median(tmpvec, vii);
    for (j=0; j<nc; j++) {
      ii = i + j*nr;
      if (idxMat[ii]) x_k[ii] = med;
    }
  }

  for (i=0; i<len; i++) {
    if (x_k[i] < minvalue) x_k[i] = minvalue;
  }

}


static void get_rowsums(x, nr, nc, ret)
double *x, *ret; /* x is vector of  stacked columns */
int nr, nc;
{
  int i, j, k=0;

  for (i=0; i<nr; i++) ret[i] = 0.0;

  for (j=0; j<nc; j++) {
    for (i=0; i<nr; i++) {
      ret[i] += x[k];
      k++; 
    }
  }

}

static void get_colsums(x, nr, nc, ret)
double *x, *ret; /* x is vector of  stacked columns */
int nr, nc;
{
  int i, j, k=0;
  double sum;

  for (j=0; j<nc; j++) {
    sum = 0.0;
    for (i=0; i<nr; i++) {
      sum += x[k];
      k++;
    }
    ret[j] = sum;
  }

}

static void transpose(v, nr, nc, ret)
double *v, *ret;
int nr, nc;
{
  int i, j, k=0;

  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      ret[k] = v[j*nr + i]; 
      k++;
    }
  }

}

static void ECM_update_h(h1, h0, w0, x_k, hnr, hnc, wnr, wnc, xnr, xnc, vM, vN, vN2)
double *h1, *h0, *w0, *x_k, *vM, *vN, *vN2;
int hnr, hnc, wnr, wnc, xnr, xnc;
{
  int i, j, k, xlen, hlen;

  xlen = xnr*xnc;
  hlen = hnr*hnc;

  /* row sums of t(w0) = col sums of w0 (t(w0) %*% mat1) */
  get_colsums(w0, wnr, wnc, vM);

  /* compute: input.k/(W0 %*% H0) */
  matrixMult(w0, wnr, wnc, h0, hnc, vN);
  for (i=0; i<xlen; i++) vN[i] = x_k[i]/vN[i];
   
  /* transpose w0 */
  transpose(w0, wnr, wnc, vN2);

  /* compute: t(w0) %*% (input.k/(W0 %*% H0)) */
  matrixMult(vN2, wnc, wnr, vN, xnc, h1);
  

  /* compute: H0*(t(W0) %*% (input.k/(W0 %*% H0)))/(t(W0) %*% mat1)
     wnc = hnr */
  k = 0;
  for (i=0; i<hnc; i++) {
    for (j=0; j<hnr; j++) {
      h1[k] = h0[k]*h1[k]/vM[j];
      k++;
    }
  }

}

static void ECM_update_w(w1, h1, w0, x_k, hnr, hnc, wnr, wnc, xnr, xnc, vM, vN, vN2)
double *h1, *w1, *w0, *x_k, *vM, *vN, *vN2;
int hnr, hnc, wnr, wnc, xnr, xnc;
{
  int i, j, k, xlen, hlen;
  double tmp;

  xlen = xnr*xnc;
  hlen = hnr*hnc;

  /* transpose h1 */
  transpose(h1, hnr, hnc, vN2);

  /* col sums of t(h1), (mat1 %*% t(h1)) */
  get_colsums(vN2, hnc, hnr, vM);

  /* compute: input.k/(W0 %*% H1) */
  matrixMult(w0, wnr, wnc, h1, hnc, vN);
  for (i=0; i<xlen; i++) vN[i] = x_k[i]/vN[i];
   
  /* compute:  (input.k/(W0 %*% H0) %*% t(h1)) */
  matrixMult(vN, xnr, xnc, vN2, hnr, w1);
  
  /* compute: W0*((input.k/(W0 %*% H1)) %*% t(H1))/(mat1 %*% t(H1)) */
  k = 0;
  for (i=0; i<wnc; i++) {
    tmp = vM[i];
    for (j=0; j<wnr; j++) {
      w1[k] = w0[k]*w1[k]/tmp;
      k++;
    }
  }

}

static int ECM_alg(str)
NMFSTR *str;
{
  /* x_k, x_k_hat, delta_denom, EM_iter get updated in str */

  int conv=0, iter, len, i, xnr, xnc, wnr, wnc, hnr, hnc, wlen, hlen;
  double *x, *logx, *x_k, *x_k_hat, minval, delta_denom, lik1, lik0; 
  double train_adj, tmp, eps, *w, *h, *w0, *h0, *vM, *vN, *vN2, *p1, *p2;
  binary *idxMat, *pb;

  x         = str->x;
  logx      = str->logx;
  x_k       = str->x_k;
  x_k_hat   = str->x_k_hat;
  w         = str->w;
  w0        = str->w0;
  h         = str->h;
  h0        = str->h0;
  xnr       = str->xnr;
  xnc       = str->xnc;
  wnr       = str->wnr;
  wnc       = str->wnc;
  hnr       = str->hnr;
  hnc       = str->hnc;
  wlen      = str->wlen;
  hlen      = str->hlen;
  vM        = str->tmpVecM;
  vN        = str->tmpVecN;
  vN2       = str->tmpVecN2;
  len       = str->xlen;
  minval    = str->EM_minvalue;
  idxMat    = str->idxMat;
  train_adj = str->train_adj;
  lik0      = 1.0; 
  eps       = str->EM_eps;
  
  copy_dVec(w, wlen, w0);
  copy_dVec(h, hlen, h0);

  for (iter=1; iter<=str->EM_maxiter; iter++) {
    replaceZeroVal_dVec(x_k, len, minval);

    /* Update W and H */
    ECM_update_h(h, h0, w0, x_k, hnr, hnc, wnr, wnc, xnr, xnc, vM, vN, vN2);  
    ECM_update_w(w, h,  w0, x_k, hnr, hnc, wnr, wnc, xnr, xnc, vM, vN, vN2);

    /* x_k_hat = W %*% H */
    matrixMult(w, wnr, wnc, h, hnc, x_k_hat);

    delta_denom = get_delta_denom(x_k_hat, len);
    lik1 = loglike(x, logx, x_k_hat, idxMat, delta_denom, len, 0);
    lik1 = sqrt(lik1/train_adj);
    tmp  = fabs((lik0 - lik1)/lik0);
    if ((iter > 1) && (tmp < eps)) {
       conv = 1;
       break;
    } 
    lik0 = lik1;
    copy_dVec(w, wlen, w0);
    copy_dVec(h, hlen, h0);

    for (i=0, pb=idxMat, p1=x_k, p2=x_k_hat; i<len; i++, pb++, p1++, p2++) {
      if (*pb) *p1 = *p2;
    }
  }

  str->delta_denom = delta_denom;
  str->EM_iter     = iter;
  str->loglike     = lik1;

  return(conv);

}

static void init_wh_rank(str, rank)
NMFSTR *str;
int rank;
{
  int n;

  str->wnc     = rank;
  str->hnr     = rank;
  str->wlen    = (str->wnr)*(str->wnc);
  str->hlen    = (str->hnr)*(str->hnc);

  if (str->w) free(str->w);
  if (str->h) free(str->h);
  if (str->w0) free(str->w0);
  if (str->h0) free(str->h0);

  str->rank    = rank;
  n            = str->wlen;
  str->w       = dVec_alloc(n, 0, 0.0);
  str->w0      = dVec_alloc(n, 0, 0.0);
  n            = str->hlen;
  str->h       = dVec_alloc(n, 0, 0.0);
  str->h0      = dVec_alloc(n, 0, 0.0);

}

static void suitor_main(str)
NMFSTR *str;
{
  /* The seed has been set i the R code before this call.
     kVec and rankVec should be ordered by kVec
  */
  int i, *kVec, *rankVec, k, rank, print, k0, xnr, xnc, kfold, xlen, conv, sumIdx;
  double minvalue, *x, *x_k, *tmpvec, train_err, test_err, *ret_train, *ret_test;
  double delta_denom, *x_k_hat, *logx, *x_k0;
  binary *idxMat;

  k0        = -1;
  x         = str->x;
  x_k       = str->x_k;
  x_k_hat   = str->x_k_hat;
  x_k0      = str->x_k0;
  logx      = str->logx;
  idxMat    = str->idxMat;
  kfold     = str->kfold;
  kVec      = str->kVec;
  rankVec   = str->rankVec;
  print     = str->print;
  xlen      = str->xlen;
  ret_train = str->train_errors;
  ret_test  = str->test_errors;
  xnr       = str->xnr;
  xnc       = str->xnc;
  tmpvec    = str->tmpVecN;
  minvalue  = str->EM_minvalue;


  for (i=0; i<str->nk; i++) {
    k    = kVec[i];
    rank = rankVec[i];
    if (print) Rprintf("rank = %d, k = %d\n", rank, k); 

    /* init objects that depend on the rank */
    init_wh_rank(str, rank);
 
    if (k != k0) {
      sumIdx = get_idxMat0(idxMat, xnr, xnc, k, kfold); 

      /* Initialize x_k */
      get_init_x_k(x, x_k, idxMat, xnr, xnc, tmpvec, minvalue);
      copy_dVec(x_k, xlen, x_k0);
    } else {
      copy_dVec(x_k0, xlen, x_k);
    }
    k0 = k;

    /* Call nmf function, Randomly assign w and h inside */
    call_nmf(x_k, str);

    /* x_k_hat is W%*%H */
    matrixMult(str->w, str->wnr, str->wnc, str->h, str->hnc, x_k_hat);

    delta_denom    = get_delta_denom(x_k_hat, xlen);
    str->train_adj = xlen - sumIdx;
    update_x_k(x_k, x_k_hat, idxMat, xlen); /* Change this */

    /* Call ECM */
    conv = ECM_alg(str);

    if (conv) {
      if (print > 2) Rprintf("EM algorithm converged in %d iterations\n", str->EM_iter);
      delta_denom = str->delta_denom;
      test_err    = loglike(x, logx, x_k_hat, idxMat, delta_denom, xlen, 1);
      train_err   = loglike(x, logx, x_k_hat, idxMat, delta_denom, xlen, 0);
    } else {
      if (print) Rprintf("EM algorithm did not converge\n");
      test_err    = DOUBLE_MISS;
      train_err   = DOUBLE_MISS;
    }

    ret_train[i] = train_err;
    ret_test[i]  = test_err;
  }

}

/* log(x) is set to 0 if x = 0 */
static void init_logx(x, logx, n)
double *x, *logx;
int n;
{
  int i;
  double tmp;

  for (i=0; i<n; i++) {
    tmp = x[i];
    if (tmp > NUMERICZERO) {
      logx[i] = log(tmp);
    } else {
      logx[i] = 0.0;
    }
  }
}


static void NMFSTR_init(str, iargs, dargs, xmat) 
NMFSTR *str;
int *iargs;
double *dargs, *xmat;
{
  int n, i;

  str->x            = NULL;
  str->x_k          = NULL;
  str->x_k_hat      = NULL;
  str->w            = NULL;
  str->w0           = NULL;
  str->h            = NULL;
  str->h0           = NULL;
  str->logx         = NULL;
  str->rankVec      = NULL;
  str->kVec         = NULL;
  str->runifVec     = NULL;
  str->cons         = NULL;
  str->h_max_index  = NULL;
  str->train_errors = NULL;
  str->test_errors  = NULL;
  str->tmpVecM      = NULL;
  str->tmpVecN      = NULL;
  str->tmpVecN2     = NULL;
  str->idxMat       = NULL;
  str->idxInt       = NULL;

  str->x           = xmat;
  str->xnr         = iargs[IARG_NR];
  str->xnc         = iargs[IARG_NC];
  str->xlen        = (str->xnr)*(str->xnc);
  str->rank        = iargs[IARG_RANK];
  str->kfold       = iargs[IARG_KFOLD];
  str->EM_maxiter  = iargs[IARG_EM_MAXITER];
  str->maxRank     = iargs[IARG_MAX_RANK];
  str->print       = iargs[IARG_PRINT];
  str->nk          = iargs[IARG_NUM_K];
  
  str->EM_eps      = dargs[DARG_EM_EPS];
  str->update_eps  = dargs[DARG_UPDATE_EPS];
  str->EM_minvalue = dargs[DARG_EM_MINVALUE];

  str->EM_iter     = 0;
  str->delta_denom = 0.0;
  str->train_adj   = 0.0;
  str->loglike     = 0.0;
  
  /* These remain fixed, others depend on rank, see init_wh_rank */
  str->wnr         = str->xnr;
  str->hnc         = str->xnc;

  str->cons_inc    = 0;
  str->cons_len    = (str->hnc)*(str->hnc - 1)/2;
  str->cons        = cVec_alloc(str->cons_len, 1, 0);
  str->h_max_index = iVec_alloc(str->hnc, 0, 0);

  /* Vector to store uniform random numbers for W and H */
  n = (str->xnr)*(str->maxRank) + (str->maxRank)*(str->xnc);
  str->runifVec = dVec_alloc(n, 0, 0.0);
  get_runif_vec(str->runifVec, n, 0.0, 1.0);

  /* Only if suitor_main will be called */
  if (str->nk) {

    n = str->xlen;
    str->idxMat = cVec_alloc(n, 0, 0);

    /* idxInt must have length at least ncols*ceil(nrows/kfold) */
    n = str->xnc*ceil((str->xnr)/(str->kfold));
    str->idxInt   = iVec_alloc(n, 0, 0);

    n             = str->xlen;
    str->logx     = dVec_alloc(n, 0, 0.0);
    init_logx(str->x, str->logx, n);
    str->x_k      = dVec_alloc(n, 0, 0.0);
    str->x_k_hat  = dVec_alloc(n, 0, 0.0);
    str->x_k0     = dVec_alloc(n, 0, 0.0);
    str->tmpVecN  = dVec_alloc(n, 0, 0.0);
    str->tmpVecN2 = dVec_alloc(n, 0, 0.0);
    n             = MAX(str->xnr, str->xnc);
    str->tmpVecM  = dVec_alloc(n, 0, 0.0);
    
  } else {
    /* init objects that depend on the rank */
    init_wh_rank(str, str->rank);
  }
}

static void NMFSTR_free(str, EMflag)
NMFSTR *str;
int EMflag;
{

  free(str->cons);
  free(str->h_max_index);
  free(str->runifVec);
  
  if (EMflag) {
    free(str->x_k);
    free(str->x_k_hat);
    free(str->x_k0);
    free(str->logx);
    free(str->w);
    free(str->h);
    free(str->w0);
    free(str->h0);
    free(str->tmpVecN);
    free(str->tmpVecN2);
    free(str->tmpVecM);

    free(str->idxMat);
    free(str->idxInt);
  }

}

void C_call_nmf(mat, iargs, dargs, retW, retH)
double *mat, *dargs, *retW, *retH;
int *iargs;
{
  NMFSTR str;

  /* For random number generation */
  GetRNGstate();

  NMFSTR_init(&str, iargs, dargs, mat);
  str.w = retW;
  str.h = retH;

  call_nmf(mat, &str);

  NMFSTR_free(&str, 0);

  PutRNGstate(); 

  return;
}

void C_call_suitor(mat, iargs, dargs, rankVec, kVec, ret_train, ret_test)
double *mat, *dargs, *ret_train, *ret_test;
int *iargs, *rankVec, *kVec;
{
  NMFSTR str;

  /* For random number generation */
  GetRNGstate();

  NMFSTR_init(&str, iargs, dargs, mat);
  str.rankVec      = rankVec;
  str.kVec         = kVec;
  str.train_errors = ret_train;
  str.test_errors  = ret_test;

  suitor_main(&str);

  NMFSTR_free(&str, 1);

  PutRNGstate(); 

  return;
}





void R_init_SUITOR(DllInfo *dll)
{
    R_registerRoutines(dll, callMethods, NULL, NULL, NULL);
}

