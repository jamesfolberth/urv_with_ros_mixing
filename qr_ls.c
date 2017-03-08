// Compile with (at least for R2015b and R2016a on GNU/Linux)
//    mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99" qr_ls.c -lmwlapack

#include <assert.h>
#include <string.h>

#include "mex.h"
#include "lapack.h"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define MINNORM_TRANSPOSE 1 // a bit faster for large problems
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

void dqr_ls(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
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
      mexErrMsgTxt("qr_ls.c: dimension 1 of A does not match dimension 1 of b.");
   
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
   tau = mxMalloc(MIN(m,n)*sizeof(*Aio));
   ptrdiff_t wrk_size = -1; // -1 gives workspace query
   wrk = mxMalloc(sizeof(*Aio));
   ptrdiff_t info;
   ptrdiff_t k = MIN(m,n);

   if ( m >= n ) { // overdetermined or square
      // Make a copy of A (used for work) (LAPACK overwrites) 
      Aio_mx = mxDuplicateArray(prhs[0]); Aio = mxGetPr(Aio_mx);

      // Do all the workspace size queries up front so we allocate wrk once
      ptrdiff_t wrks_qrf=-1, wrks_mqr=-1;
      dgeqrf(&m,&n,Aio,&m,tau,wrk,&wrks_qrf,&info);
      wrks_qrf = (ptrdiff_t)wrk[0];
      dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrk_size,&info);
      wrks_mqr = (ptrdiff_t)wrk[0];
      wrk_size = MAX(wrks_qrf, wrks_mqr);
      mxFree(wrk);
      wrk = mxMalloc(wrk_size*sizeof(*Aio));
#if DEBUG > 0
      printf("wrk_size (DGEQRF) = %d\n", wrks_qrf);
      printf("wrk_size (DORMQR) = %d\n", wrks_mqr);
#endif
      
      // Perform QR factorization A=QR
      dgeqrf(&m,&n,Aio,&m,tau,wrk,&wrks_qrf,&info);
      if ( info < 0 )
         mexErrMsgTxt("qr_ls.c: DGEQRF failed.");

      // Apply Q^T to b
      //TODO This gives the right answer, but applies the full Q^T;
      //     we want just the first MIN(m,n) columns of Q to be applied in Q^T
      //     Either way, this is still considerably faster than building the
      //     trimmed Q and DGEMM'ing it
      dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrks_mqr,&info);
      if ( info < 0 )
         mexErrMsgTxt("qr_ls.c: DORMQR failed.");

      // Solve triangular system Rx=(Q^Tb)
      dtrtrs("U","N","N",&n,&n_rhs,Aio,&m,bio,&m,&info);
      if ( info != 0 )
         mexErrMsgTxt("qr_ls.c: DTRTRS failed.  Is A rank-deficient?");

      // copy "x" from b to x
      for (ptrdiff_t j=0; j<n_rhs; ++j) {
         //for (i=0; i<n; ++i) {
         //   *(x+i+n*j) = *(bio+i+m*j);
         //}
         memcpy(x+n*j,bio+m*j,n*sizeof(*bio));
      }

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
            mexErrMsgTxt("qr_ls.c: mode should be a string ('basic' or 'minnorm')");

         mode = mxArrayToString(prhs[2]);
         if ( mode == NULL ) mexErrMsgTxt("qr_ls.c: error converting mode to C string");
         if ( strlen(mode) < 5 ) mexErrMsgTxt("qr_ls.c: mode should be either "
               "'basic' or 'minnorm'.");
      }
      
      // if basic solution requested, do QR on trimmed A (A(1:m,1:m)), solve, and
      // finally pad solution with zeros.
      if ( strcmp(mode, "basic") == 0 ) {
         // Make a copy of A and copy b into x (used for work) (LAPACK overwrites) 
         // TODO: for the 'basic' solution, we only need entries 1:m,1:m instead of 1:m,1:n
         Aio_mx = mxDuplicateArray(prhs[0]); Aio = mxGetPr(Aio_mx);

         // Do all the workspace size queries up front so we allocate wrk once
         ptrdiff_t wrks_qrf=-1, wrks_mqr=-1;
         dgeqrf(&m,&m,Aio,&m,tau,wrk,&wrks_qrf,&info);
         wrks_qrf = (ptrdiff_t)wrk[0];
         dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrk_size,&info);
         wrks_mqr = (ptrdiff_t)wrk[0];
         wrk_size = MAX(wrks_qrf, wrks_mqr);
         mxFree(wrk);
         wrk = mxMalloc(wrk_size*sizeof(*Aio));
#if DEBUG > 0
         printf("wrk_size (DGEQRF) = %d\n", wrks_qrf);
         printf("wrk_size (DORMQR) = %d\n", wrks_mqr);
#endif

         // Perform QR factorization A=QR
         dgeqrf(&m,&m,Aio,&m,tau,wrk,&wrks_qrf,&info);
         if ( info < 0 )
            mexErrMsgTxt("qr_ls.c: DGEQRF failed.");

         // Apply Q^T to b
         dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,bio,&m,wrk,&wrks_mqr,&info);
         if ( info < 0 )
            mexErrMsgTxt("qr_ls.c: DORMQR failed.");

         // Solve triangular system Rx=(Q^Tb)
         dtrtrs("U","N","N",&m,&n_rhs,Aio,&m,bio,&m,&info);
         if ( info != 0 )
            mexErrMsgTxt("qr_ls.c: DTRTRS failed.  Is A rank-deficient?");

         // copy "x" from b to x and zero pad
         for (ptrdiff_t j=0; j<n_rhs; ++j) {
            memcpy(x+n*j,bio+m*j,m*sizeof(*bio));
            memset(x+n*j+m,0,(n-m)*sizeof(*x));
         }
      }

#if MINNORM_TRANSPOSE == 1
      // if minimum norm solution requested, do QR on A^T.  Q^T will be the same size
      // as A, so this is going to be slower than the basic solution.
      else if ( strcmp(mode, "minnorm") == 0 ) {
         // Make a copy of A^T and copy b into x (used for work) (LAPACK overwrites) 
         mexCallMATLAB(1, &Aio_mx, 1, (mxArray **) &prhs[0], "transpose");
         Aio = mxGetPr(Aio_mx); // Note: Aio is actually A^T

         // Do all the workspace size queries up front so we allocate wrk once
         ptrdiff_t wrks_qrf=-1, wrks_mqr=-1;
         dgeqrf(&n,&m,Aio,&n,tau,wrk,&wrks_qrf,&info);
         wrks_qrf = (ptrdiff_t)wrk[0];
         dormqr("L","T",&m,&n_rhs,&k,Aio,&m,tau,x,&n,wrk,&wrk_size,&info);
         wrks_mqr = (ptrdiff_t)wrk[0];
         wrk_size = MAX(wrks_qrf, wrks_mqr);
         mxFree(wrk);
         wrk = mxMalloc(wrk_size*sizeof(*Aio));
#if DEBUG > 0
         printf("wrk_size (DGEQRF) = %d\n", wrks_qrf);
         printf("wrk_size (DORMQR) = %d\n", wrks_mqr);
#endif
         
         // Perform QR factorization A^T=QR
         dgeqrf(&n,&m,Aio,&n,tau,wrk,&wrks_qrf,&info);
         if ( info < 0 )
            mexErrMsgTxt("qr_ls.c: DGEQRF failed.");
         
         // Solve triangular system R^T(Q^Tx)=R^Ty=b
         dtrtrs("U","T","N",&m,&n_rhs,Aio,&n,bio,&m,&info);
         if ( info != 0 )
            mexErrMsgTxt("qr_ls.c: DTRTRS failed.  Is A rank-deficient?");

         // Copy components of b into x and zero uninitialized section
         for (ptrdiff_t j=0; j<n_rhs; ++j) {
            memcpy(x+n*j,bio+m*j,m*sizeof(*bio));
            memset(x+n*j+m,0,(n-m)*sizeof(*x));
         }

         // Find x via x = Qy 
         dormqr("L","N",&n,&n_rhs,&k,Aio,&n,tau,x,&n,wrk,&wrk_size,&info);
         if ( info < 0 )
            mexErrMsgTxt("qr_ls.c: DORMQR failed.");

      }
#else // MINNORM_TRANSPOSE != 1
      // if minimum norm solution requested, do LQ on A.  Q will be the same size
      // as A, so this is going to be slower than the basic solution.
      else if ( strcmp(mode, "minnorm") == 0 ) {
         // Make a copy of A and copy b into x (used for work) (LAPACK overwrites) 
         // TODO: for the 'basic' solution, we only need entries 1:m,1:m instead of 1:m,1:n
         Aio_mx = mxDuplicateArray(prhs[0]); Aio = mxGetPr(Aio_mx);

         // Do all the workspace size queries up front so we allocate wrk once
         ptrdiff_t wrks_lqf=-1, wrks_mlq=-1;
         dgelqf(&m,&n,Aio,&m,tau,wrk,&wrks_lqf,&info);
         wrks_lqf = (ptrdiff_t)wrk[0];
         dormlq("L","T",&n,&n_rhs,&k,Aio,&m,tau,x,&n,wrk,&wrk_size,&info);
         wrks_mlq = (ptrdiff_t)wrk[0];
         wrk_size = MAX(wrks_lqf, wrks_mlq);
         mxFree(wrk);
         wrk = mxMalloc(wrk_size*sizeof(*Aio));
#if DEBUG > 0
         printf("wrk_size (DGELQF) = %d\n", wrks_lqf);
         printf("wrk_size (DORMLQ) = %d\n", wrks_mlq);
#endif
         
         // Perform LQ factorization A=LQ
         dgelqf(&m,&n,Aio,&m,tau,wrk,&wrks_lqf,&info);
         if ( info < 0 )
            mexErrMsgTxt("qr_ls.c: DGELQF failed.");
         
         
         // Solve triangular system L(Qx)=Ly=b
         dtrtrs("L","N","N",&m,&n_rhs,Aio,&m,bio,&m,&info);
         if ( info != 0 )
            mexErrMsgTxt("qr_ls.c: DTRTRS failed.  Is A rank-deficient?");
         
         // Copy components of b into x and zero uninitialized section
         for (ptrdiff_t j=0; j<n_rhs; ++j) {
            memcpy(x+n*j,bio+m*j,m*sizeof(*bio));
            memset(x+n*j+m,0,(n-m)*sizeof(*x));
         }

         // Find x via x = Q^Ty 
         dormlq("L","T",&n,&n_rhs,&k,Aio,&m,tau,x,&n,wrk,&wrk_size,&info);
         if ( info < 0 )
            mexErrMsgTxt("qr_ls.c: DORMLQ failed.");

      }
#endif

      else 
         mexErrMsgTxt("qr_ls.c: mode should be either 'basic' or 'minnorm'.");

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
      mexErrMsgTxt("qr_ls.c: qr_ls expects 2 or 3 input arguments.");

   if ( nlhs > 1 )
      mexWarnMsgTxt("qr_ls.c: qr_ls uses only 1 output argument.");
   
   bool isd = mxIsDouble(prhs[0]);
   bool iss = mxIsSingle(prhs[0]);
   if ( !(isd || iss) || mxIsComplex(prhs[0]) )
      mexErrMsgTxt("qr_ls.c: qr_ls only works on dense, real-valued matrices.");
   else if ( isd )
      dqr_ls(nlhs, plhs, nrhs, prhs);
   else if ( iss )
      //sqr_ls(nlhs, plhs, nrhs, prhs);
      mexErrMsgTxt("qr_ls.c: single not supported.");

   return;
}
