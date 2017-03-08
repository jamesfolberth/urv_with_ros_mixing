// Compile with (at least for R2015b and R2016a on GNU/Linux)
//    mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99" qrcp_ls.c -lmwlapack

#include <assert.h>
#include <string.h>

#include "mex.h"
#include "lapack.h"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

//#define MINNORM_TRANSPOSE 1 // LAPACK doesn't have a pivoted LQ, 
                              // so we always transpose and use _geqp3
#define DEBUG 0

#if DEBUG > 0
void print_array(mxArray *A) {
   mexCallMATLAB(0,NULL,1,&A,"disp");
}

void print_array_label(mxArray *A, const char *label) {
   assert(label[strlen(label)+1] == '\0');
   char buf[strlen(label)+3];
   strncpy(buf, label, strlen(label)+1);
   strcat(buf, ":\n");
   printf(buf);
   mexCallMATLAB(0,NULL,1,&A,"disp");
}
#endif


void dqrcp_ls(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   double *x,*A,*b;
   double *Aio,*bio; // in/out for LAPACK
   mxArray *Aio_mx,*bio_mx; // in/out for LAPACK
   double *tau, *wrk; // used for LAPACK
   ptrdiff_t dims[2];
   ptrdiff_t m,n;
   ptrdiff_t n_rhs; // number of RHS vectors
   
   m = mxGetM(prhs[0]);
   n = mxGetN(prhs[0]);
   n_rhs = mxGetN(prhs[1]);
   if ( m != mxGetM(prhs[1]) )
      mexErrMsgTxt("qrcp_ls.c: dimension 1 of A does not match dimension 1 of b.");
   
   mxClassID class = mxGetClassID(prhs[0]);
   assert( class == mxDOUBLE_CLASS || mxSINGLE_CLASS );

   A = mxGetPr(prhs[0]);
   b = mxGetPr(prhs[1]);
   // Copy b into x (used for work) (LAPACK overwrites) 
   //Aio_mx = mxDuplicateArray(prhs[0]); Aio = mxGetPr(Aio_mx);
   bio_mx = mxDuplicateArray(prhs[1]); bio = mxGetPr(bio_mx);
   
   dims[0] = n; dims[1] = n_rhs;
   plhs[0] = mxCreateUninitNumericArray(2, dims, class, mxREAL); // makes output vector(s)
   x = mxGetPr(plhs[0]);
   
   // variables for LAPACK work
   ptrdiff_t *pvt, *ipvt;
   tau = mxMalloc(MIN(m,n)*sizeof(*Aio));
   ptrdiff_t wrk_size = -1; // -1 gives workspace query
   wrk = mxMalloc(sizeof(*Aio));
   ptrdiff_t info;
   ptrdiff_t k = MIN(m,n);

   if ( m >= n ) { // overdetermined or square
      // Make a copy of A (used for work) (LAPACK overwrites) 
      Aio_mx = mxDuplicateArray(prhs[0]); Aio = mxGetPr(Aio_mx);
      
      // allocate pivot and inverse pivot vectors
      pvt = mxCalloc(n,sizeof(*pvt)); 
      ipvt = mxMalloc(n*sizeof(*ipvt));

      // Do all the workspace size queries up front so we allocate wrk once
      ptrdiff_t wrks_qp3=-1, wrks_mqr=-1;
      dgeqp3(&m,&n,Aio,&m,pvt,tau,wrk,&wrks_qp3,&info);
      wrks_qp3 = (ptrdiff_t)wrk[0];
      dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrk_size,&info);
      wrks_mqr = (ptrdiff_t)wrk[0];
      wrk_size = MAX(wrks_qp3, wrks_mqr);
      mxFree(wrk);
      wrk = mxMalloc(wrk_size*sizeof(*Aio));
#if DEBUG > 0
      printf("wrk_size (DGEQP3) = %d\n", wrks_qp3);
      printf("wrk_size (DORMQR) = %d\n", wrks_mqr);
#endif
      
      // Perform pivoted QR factorization A(:,pvt)=QR
      dgeqp3(&m,&n,Aio,&m,pvt,tau,wrk,&wrks_qp3,&info);
      if ( info < 0 )
         mexErrMsgTxt("qrcp_ls.c: DGEQP3 failed.");

      // Switch to 0-based indexing and populate inverse pivot vector
      for (ptrdiff_t i=0; i<n; ++i) {
         pvt[i] -= 1; // pvt is 1-based coming out of LAPACK
         ipvt[pvt[i]] = i;
      }
      
      // Apply Q^T to b
      //TODO This gives the right answer, but applies the full Q^T;
      //     we want just the first MIN(m,n) columns of Q to be applied in Q^T
      //     Either way, this is still considerably faster than building the
      //     trimmed Q and DGEMM'ing it
      dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrks_mqr,&info);
      if ( info < 0 )
         mexErrMsgTxt("qrcp_ls.c: DORMQR failed.");

      // Solve triangular system Rx=(Q^Tb)
      dtrtrs("U","N","N",&n,&n_rhs,Aio,&m,bio,&m,&info);
      if ( info != 0 )
         mexErrMsgTxt("qrcp_ls.c: DTRTRS failed.  Is A rank-deficient?");
      
      // copy "x" from b to x and inverse pivot
      for (ptrdiff_t j=0; j<n_rhs; ++j) {
         for (ptrdiff_t i=0; i<n; ++i) {
            *(x+i+n*j) = *(bio+ipvt[i]+m*j);
         }
      }
      
      mxFree(pvt);
      mxFree(ipvt);
   } 

   else { // underdetermined
      // if present, get mode == 'basic' or 'minnorm'
      char *mode;
      if ( nrhs == 2 ) {
         mode = mxMalloc(6*sizeof(char));
         strcpy(mode, "basic");
      }
      else {
         if ( !mxIsChar(prhs[2]) )
            mexErrMsgTxt("qrcp_ls.c: mode should be a string ('basic' or 'minnorm')");

         mode = mxArrayToString(prhs[2]);
         if ( mode == NULL ) mexErrMsgTxt("qrcp_ls.c: error converting mode to C string");
         if ( strlen(mode) < 5 ) mexErrMsgTxt("qrcp_ls.c: mode should be either "
               "'basic' or 'minnorm'.");
      }
      
      // if basic solution requested, do QR on trimmed A (A(1:m,1:m)), solve, and
      // finally pad solution with zeros.
      if ( strcmp(mode, "basic") == 0 ) {
         // Make a copy of A and copy b into x (used for work) (LAPACK overwrites) 
         // TODO: for the 'basic' solution, we only need entries 1:m,1:m instead of 1:m,1:n
         Aio_mx = mxDuplicateArray(prhs[0]); Aio = mxGetPr(Aio_mx);

         // allocate pivot and inverse pivot vectors
         pvt = mxCalloc(m,sizeof(*pvt)); 
         ipvt = mxMalloc(m*sizeof(*ipvt));

         // Do all the workspace size queries up front so we allocate wrk once
         ptrdiff_t wrks_qp3=-1, wrks_mqr=-1;
         dgeqp3(&m,&m,Aio,&m,pvt,tau,wrk,&wrks_qp3,&info);
         wrks_qp3 = (ptrdiff_t)wrk[0];
         dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrk_size,&info);
         wrks_mqr = (ptrdiff_t)wrk[0];
         wrk_size = MAX(wrks_qp3, wrks_mqr);
         mxFree(wrk);
         wrk = mxMalloc(wrk_size*sizeof(*Aio));
#if DEBUG > 0
         printf("wrk_size (DGEQP3) = %d\n", wrks_qp3);
         printf("wrk_size (DORMQR) = %d\n", wrks_mqr);
#endif

         // Perform QR factorization A(1:m,pvt(1:m))=QR
         dgeqp3(&m,&m,Aio,&m,pvt,tau,wrk,&wrks_qp3,&info);
         if ( info < 0 )
            mexErrMsgTxt("qrcp_ls.c: DGEQP3 failed.");

         // Switch to 0-based indexing and populate inverse pivot vector
         for (ptrdiff_t i=0; i<m; ++i) {
            pvt[i] -= 1; // pvt is 1-based coming out of LAPACK
            ipvt[pvt[i]] = i;
         }

         // Apply Q^T to b
         dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrks_mqr,&info);
         if ( info < 0 )
            mexErrMsgTxt("qrcp_ls.c: DORMQR failed.");

         // Solve triangular system Rx=(Q^Tb)
         dtrtrs("U","N","N",&m,&n_rhs,Aio,&m,bio,&m,&info);
         if ( info != 0 )
            mexErrMsgTxt("qrcp_ls.c: DTRTRS failed.  Is A rank-deficient?");

         // copy "x" from b to x and zero pad
         for (ptrdiff_t j=0; j<n_rhs; ++j) {
            for (ptrdiff_t i=0; i<m; ++i) {
               *(x+i+n*j) = *(bio+ipvt[i]+m*j);
            }
            memset(x+n*j+m,0,(n-m)*sizeof(*x));
         }

         mxFree(pvt);
         mxFree(ipvt);
      }

      // if basicfull solution requested, do QRCP on fullr A, solve trimmed system, and
      // and place zeros where appropriate
      else if ( strcmp(mode, "basicfull") == 0 ) {
         Aio_mx = mxDuplicateArray(prhs[0]); Aio = mxGetPr(Aio_mx);

         // allocate pivot and inverse pivot vectors
         pvt = mxCalloc(n,sizeof(*pvt)); 
         ipvt = mxMalloc(n*sizeof(*ipvt));

         // Do all the workspace size queries up front so we allocate wrk once
         ptrdiff_t wrks_qp3=-1, wrks_mqr=-1;
         dgeqp3(&m,&n,Aio,&m,pvt,tau,wrk,&wrks_qp3,&info);
         wrks_qp3 = (ptrdiff_t)wrk[0];
         dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrk_size,&info);
         wrks_mqr = (ptrdiff_t)wrk[0];
         wrk_size = MAX(wrks_qp3, wrks_mqr);
         mxFree(wrk);
         wrk = mxMalloc(wrk_size*sizeof(*Aio));
#if DEBUG > 0
         printf("wrk_size (DGEQP3) = %d\n", wrks_qp3);
         printf("wrk_size (DORMQR) = %d\n", wrks_mqr);
#endif

         // Perform QR factorization A(:,pvt(1:n))=QR
         dgeqp3(&m,&n,Aio,&m,pvt,tau,wrk,&wrks_qp3,&info);
         if ( info < 0 )
            mexErrMsgTxt("qrcp_ls.c: DGEQP3 failed.");

         // Switch to 0-based indexing and populate inverse pivot vector
         for (ptrdiff_t i=0; i<n; ++i) {
            pvt[i] -= 1; // pvt is 1-based coming out of LAPACK
            ipvt[pvt[i]] = i;
         }

         // Apply Q^T to b
         dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrks_mqr,&info);
         if ( info < 0 )
            mexErrMsgTxt("qrcp_ls.c: DORMQR failed.");

         // Solve triangular system Rx=(Q^Tb)
         dtrtrs("U","N","N",&m,&n_rhs,Aio,&m,bio,&m,&info);
         if ( info != 0 )
            mexErrMsgTxt("qrcp_ls.c: DTRTRS failed.  Is A rank-deficient?");

         // copy "x" from b to x and zero pad
         for (ptrdiff_t j=0; j<n_rhs; ++j) {
            for (ptrdiff_t i=0; i<n; ++i) {
               double val = (ipvt[i]>=m) ? 0. : *(bio+ipvt[i]+m*j);
               *(x+i+n*j) = val;
            }
         }

         mxFree(pvt);
         mxFree(ipvt);
      }

      // if minimum norm solution requested, do QR on A^T.  Q^T will be the same size
      // as A, so this is going to be slower than the basic solution.
      else if ( strcmp(mode, "minnorm") == 0 ) {
         // Make a copy of A^T and copy b into x (used for work) (LAPACK overwrites) 
         mexCallMATLAB(1, &Aio_mx, 1, (mxArray **) &prhs[0], "transpose");
         Aio = mxGetPr(Aio_mx); // Note: Aio is actually A^T

         // allocate pivot vector (inverse is not needed)
         pvt = mxCalloc(m,sizeof(*pvt)); 

         // Do all the workspace size queries up front so we allocate wrk once
         ptrdiff_t wrks_qp3=-1, wrks_mqr=-1;
         dgeqp3(&n,&m,Aio,&n,pvt,tau,wrk,&wrks_qp3,&info);
         wrks_qp3 = (ptrdiff_t)wrk[0];
         dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,x,&n,wrk,&wrk_size,&info);
         wrks_mqr = (ptrdiff_t)wrk[0];
         wrk_size = MAX(wrks_qp3, wrks_mqr);
         mxFree(wrk);
         wrk = mxMalloc(wrk_size*sizeof(*Aio));
#if DEBUG > 0
         printf("wrk_size (DGEQP3) = %d\n", wrks_qp3);
         printf("wrk_size (DORMQR) = %d\n", wrks_mqr);
#endif
         
         // Perform QR factorization A^T(:,pvt)=QR
         dgeqp3(&n,&m,Aio,&n,pvt,tau,wrk,&wrks_qp3,&info);
         if ( info < 0 )
            mexErrMsgTxt("qrcp_ls.c: DGEQP3 failed.");

         // Switch to 0-based indexing
         for (ptrdiff_t i=0; i<m; ++i) {
            pvt[i] -= 1; // pvt is 1-based coming out of LAPACK
         }
         
         // Pivot RHS b into a copy of b
         double *biopvt;
         biopvt = mxMalloc(m*n_rhs*sizeof(*biopvt));
         for (ptrdiff_t j=0; j<n_rhs; ++j) {
            for (ptrdiff_t i=0; i<m; ++i) {
               *(biopvt+i+m*j) = *(bio+pvt[i]+m*j);
            }
         }
         
         // Solve triangular system R^T(Q^Tx)=R^Ty=b
         dtrtrs("U","T","N",&m,&n_rhs,Aio,&n,biopvt,&m,&info);
         if ( info != 0 )
            mexErrMsgTxt("qrcp_ls.c: DTRTRS failed.  Is A rank-deficient?");

         // Copy components of b into x and zero uninitialized section
         for (ptrdiff_t j=0; j<n_rhs; ++j) {
            memcpy(x+n*j,biopvt+m*j,m*sizeof(*biopvt));
            memset(x+n*j+m,0,(n-m)*sizeof(*x));
         }

         // Find x via x = Qy 
         dormqr("L","N",&n,&n_rhs,&k,Aio,&n,tau,x,&n,wrk,&wrk_size,&info);
         if ( info < 0 )
            mexErrMsgTxt("qrcp_ls.c: DORMQR failed.");

         mxFree(pvt);

      }

      else 
         mexErrMsgTxt("qrcp_ls.c: mode should be either 'basic' or 'minnorm'.");

      mxFree(mode);
   }

   //print_array_label(bio_mx,"bio");
   //print_array_label(plhs[0],"x");
   
   mxFree(tau);
   mxFree(wrk);

   return;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   // Check inputs
   if ( nrhs < 2 || nrhs > 3 )
      mexErrMsgTxt("qrcp_ls.c: qrcp_ls expects 2 or 3 input arguments.");

   if ( nlhs > 1 )
      mexWarnMsgTxt("qrcp_ls.c: qrcp_ls uses only 1 output argument.");
   
   bool isd = mxIsDouble(prhs[0]);
   bool iss = mxIsSingle(prhs[0]);
   if ( !(isd || iss) || mxIsComplex(prhs[0]) )
      mexErrMsgTxt("qrcp_ls.c: qrcp_ls only works on dense, real-valued matrices.");
   else if ( isd )
      dqrcp_ls(nlhs, plhs, nrhs, prhs);
   else if ( iss )
      //sqrcp_ls(nlhs, plhs, nrhs, prhs);
      mexErrMsgTxt("qrcp_ls.c: single not supported.");

   return;
}
