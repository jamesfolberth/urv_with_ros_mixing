// Compile with (at least for R2015b and R2016a on GNU/Linux)
//    mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99" qr_fact_only.c -lmwlapack

#include <assert.h>
#include <string.h>

#include "mex.h"
#include "lapack.h"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define DEBUG 0

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   double *A = mxGetPr(prhs[0]);
   int do_pivot = (int) mxGetScalar(prhs[1]); // 0 for no pivot, 1 for pivot

   ptrdiff_t m = mxGetM(prhs[0]);
   ptrdiff_t n = mxGetN(prhs[0]);

   // variables for LAPACK work
   double *tau, *wrk;
   tau = mxMalloc(MIN(m,n)*sizeof(*A));
   ptrdiff_t wrk_size = -1; // -1 gives workspace query
   wrk = mxMalloc(sizeof(*A));
   ptrdiff_t info;
   ptrdiff_t k = MIN(m,n);
   

   if ( do_pivot == 0 ) {
      ptrdiff_t wrks_qrf=-1;
      dgeqrf(&m,&n,A,&m,tau,wrk,&wrks_qrf,&info);
      wrks_qrf = (ptrdiff_t)wrk[0];
      mxFree(wrk);
      wrk = mxMalloc(wrks_qrf*sizeof(*A));
#if DEBUG > 0
      printf("wrk_size (DGEQRF) = %d\n", wrks_qrf);
#endif
      
      // Perform QR factorization A=QR
      dgeqrf(&m,&n,A,&m,tau,wrk,&wrks_qrf,&info);
      if ( info < 0 )
         mexErrMsgTxt("qr_fact_only.c: DGEQRF failed.");

   }
   else if ( do_pivot == 1 ) {
      ptrdiff_t *pvt;
      pvt = mxCalloc(n,sizeof(*pvt)); 

      // Do all the workspace size queries up front so we allocate wrk once
      ptrdiff_t wrks_qp3=-1;
      dgeqp3(&m,&n,A,&m,pvt,tau,wrk,&wrks_qp3,&info);
      wrks_qp3 = (ptrdiff_t)wrk[0];
      mxFree(wrk);
      wrk = mxMalloc(wrks_qp3*sizeof(*A));
#if DEBUG > 0
      printf("wrk_size (DGEQP3) = %d\n", wrks_qp3);
#endif
      
      // Perform pivoted QR factorization A(:,pvt)=QR
      dgeqp3(&m,&n,A,&m,pvt,tau,wrk,&wrks_qp3,&info);
      if ( info < 0 )
         mexErrMsgTxt("qr_fact_only.c: DGEQP3 failed.");

      mxFree(pvt);
   }
   else
      mexErrMsgTxt("bogus input");
   
   mxFree(tau);
   mxFree(wrk);

   return;
}
