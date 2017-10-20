
#include <stdio.h>
#include <stdarg.h> // vargs
#include <math.h> // sqrt
#include <string.h>
#include <ctype.h> // tolower
#include <time.h> // clock, CLOCKS_PER_SEC

#include <fftw3.h> // This uses the system FFTW (check with ldd)

#include "mex.h"

//int mxUnshareArray(mxArray *, int);

#define __USE_OMP__ 0
#ifdef __USE_OMP__
#include <omp.h>
#endif

#ifndef DEBUG
#define DEBUG 0 // debug uses the same verbosity levels as opt.verbosity
#endif

#ifndef DEBUG_OR_VERB
#define DEBUG_OR_VERB(lvl) \
   (DEBUG >= lvl || opt->verbosity >= lvl)
#endif

// http://stackoverflow.com/a/1644898
#ifndef DEBUG_PRINTF 
#define DEBUG_PRINTF(fmt, ...) \
	do { fprintf(stderr, "%s:%d:%s(): " fmt, "mix_fftw.c", \
   __LINE__, __func__, __VA_ARGS__); } while (0)
#endif

// View all the stuff in libmx.so, libmex.so
// http://stackoverflow.com/a/1620583
// readelf -Ws /usr/local/MATLAB/R2015b/bin/glnxa64/libmex.so | grep "mex"
// readelf -Ws /usr/local/MATLAB/R2015b/bin/glnxa64/libmx.so | grep "mx"
// readelf -Ws /usr/local/MATLAB/R2015b/bin/glnxa64/libmx.so | awk "{print $8}" | grep -i "trans"

// options from the mex interface
typedef struct mix_opt {
   bool transpose;
   char *side;
   char *transform;
   int plan_rigor;
   char *wisdom_file;
   bool wisdom_only;
	bool use_diags;
   int verbosity;
   int nthreads;
   size_t block_size;
} mix_opt;

// use this like printf, with a format and optional arguments
// http://stackoverflow.com/a/1485819
// http://embeddedgurus.com/stack-overflow/2008/12/efficient-c-tips-5-make-local-functions-static/
static void mex_error(char *fmt, ...) {
   va_list args;
   char *msg = mxMalloc(1024*sizeof(char));

   va_start(args, fmt);
   sprintf(msg,"mix_fftw.c: "); // hard code the file basename
                                // could do something like http://stackoverflow.com/a/8488201
   vsprintf(msg+strlen("mix_fftw.c: "), fmt, args);
   va_end(args);

   mexErrMsgTxt(msg);
}
static void mex_warning(char *fmt, ...) {
   va_list args;
   char *msg = mxMalloc(1024*sizeof(char));

   va_start(args, fmt);
   sprintf(msg,"mix_fftw.c: "); // hard code the file basename
                                // could do something like http://stackoverflow.com/a/8488201
   vsprintf(msg+strlen("mix_fftw.c: "), fmt, args);
   va_end(args);

   mexWarnMsgTxt(msg);
}


void export_wisdom(mix_opt *opt) {
   if ( strcmp(opt->wisdom_file,"") ) {
      FILE *f = fopen(opt->wisdom_file, "w");
      if ( f ) {
         if (DEBUG_OR_VERB(1))
            DEBUG_PRINTF("Exporting wisdom to file %s\n", opt->wisdom_file);
         // this concatenates new wisdom to any wisdom loaded from earlier.
         // If wisdom already exists for a specific size with a more patient
         // planner, the more rigorous wisdom is saved (and used when imported).
         clock_t end, begin=clock();
         fftw_export_wisdom_to_file(f);
         int s = ferror(f);
         if ( s != 0 )
            mex_warning("Wisdom export may have failed; file error indicator is "
                  "set for %s", opt->wisdom_file);
         s = fclose(f);
         if ( s != 0 )
            mex_warning("Wisdom export sucessful, but there was an error "
                  "closing the file %s\n", opt->wisdom_file);
         end = clock();
         if (DEBUG_OR_VERB(2))
            DEBUG_PRINTF("Exporting wisdom complete.  Wall time = %es.\n", (double)(end-begin) /
                  CLOCKS_PER_SEC);
      }
      else
         // even if issues, try to continue
         mex_warning("Error exporting wisdom to file %s\n", opt->wisdom_file);
   }
}

void import_wisdom(mix_opt *opt) {
   if ( strcmp(opt->wisdom_file,"") ) {
      FILE *f = fopen(opt->wisdom_file, "r");
      if ( f ) {
         if (DEBUG_OR_VERB(1))
            DEBUG_PRINTF("Importing wisdom from file %s\n", opt->wisdom_file);
         clock_t end, begin=clock();
         //TODO: this occasionally segfaults.  Why?!
         //      Even ..._from_filename segfaults.
         //      Restarting matlab seems to "fix" it.
         int s = fftw_import_wisdom_from_file(f);
         if ( s == 0 )
            mex_warning("Error importing wisdom from file %s\n", opt->wisdom_file);
         s = fclose(f);
         if ( s != 0 )
            mex_warning("Wisdom import sucessful, but there was an error "
                  "closing the file %s\n", opt->wisdom_file);
         end = clock();
         if (DEBUG_OR_VERB(2))
            DEBUG_PRINTF("Importing wisdom complete.  Wall time = %es.\n", (double)(end-begin) /
                  CLOCKS_PER_SEC);
      }
      else
         // Assume all is good and continue.
         // If the user specifies a new filename (i.e. file doesn't yet exist),
         // we try to read the wisdom file before it exists; this shouldn't be an 
         // issue, so we'll just continue here without warning or error.
         //mex_warning("error importing wisdom from file %s.\nDoes the file exist"
         //      " and contain FFTW wisdom?\n", opt->wisdom_file);
         if (DEBUG_OR_VERB(1))
            DEBUG_PRINTF("Issue importing wisdom from file %s.  This usually occurs when"
                  " you try to import wisdom from a new, non-existant file; in this"
                  " case, import_wisdom() will fail, we'll generate wisdom, and then "
                  "we'll export_wisdom() to the file.\n",
                  opt->wisdom_file);
   }
}


void mix_cols(double *A, double *darr, size_t m, size_t n, size_t n_its, mix_opt *opt) {
  	
	unsigned it,i,j,mj; 
   fftw_plan p;
   double *darr_it;
	double nrm_1=1,nrm_mid=1,nrm_end=1.;
   bool nrm_before=false;

   // Order of F*D_i
   bool darr_before_F, darr_reversed, F_transpose;
   bool right = strncmp(opt->side,"r",1)==0;
   if ( (right && !opt->transpose) || (!right && opt->transpose) ) { // equivalent to V^T*A
      darr_before_F = false;
      darr_reversed = true;
      F_transpose = true;
   } else { // equivalent to V*A
      darr_before_F = true;
      darr_reversed = false;
      F_transpose = false;
   }
   
	// Plan a real-to-real transform
  	// The planner syntax is 
  	// fftw_plan fftw_plan_many_r2r(int rank, const int *n, int howmany,
   //                           double *in, const int *inembed,
   //                           int istride, int idist,
   //                           double *out, const int *onembed,
   //                           int ostride, int odist,
   //                           const fftw_r2r_kind *kind, unsigned flags);
	int trans_lens[1] = {m};

	fftw_r2r_kind kind[1];
   bool do_dct = strncmp(opt->transform,"dct",4)==0;
   bool do_idct = strncmp(opt->transform,"idct",5)==0;
   bool do_dct_type = do_dct || do_idct;
   if ( do_dct_type ) {
      if ( (!F_transpose && do_dct) || (F_transpose && do_idct) ) { // DCT-II
         kind[0] = FFTW_REDFT10; 
         nrm_1 = 1./(2.*sqrt(m));
         nrm_mid = 1./sqrt(2.*m);
         nrm_end = 1./sqrt(2.*m);
         nrm_before = false;
      } else { // DCT-II
         kind[0] = FFTW_REDFT01;
         nrm_1 = 1./sqrt(m);
         nrm_mid = 1./sqrt(2.*m);
         nrm_end = 1./sqrt(2.*m);
         nrm_before = true;
      }
   }
   else
      mex_error("unknown transform type: %s", opt->transform);
   
   unsigned fftw_flags=0;
   switch (opt->plan_rigor) {
      case 1:
         fftw_flags |= FFTW_ESTIMATE;
      case 2:
         fftw_flags |= FFTW_MEASURE;
      case 3:
         fftw_flags |= FFTW_PATIENT;
      case 4:
         fftw_flags |= FFTW_EXHAUSTIVE;
   }

   // Try to get the plan from wisdom
   // If using FFTW_ESTIMATE, this forms the plan if wisdom is not present; if 
   // wisdom is present for a more patient planner, then it uses the more patient plan.
   clock_t end,begin=clock();
   p = fftw_plan_many_r2r(1, trans_lens, (int)n, A, NULL, 1, (int)m, A, NULL, 1, (int)m,
      kind, fftw_flags | FFTW_WISDOM_ONLY);
   end=clock();
   if (DEBUG_OR_VERB(2))
      DEBUG_PRINTF("FFTW_WISDOM_ONLY time = %e.\n", (double)(end-begin)/CLOCKS_PER_SEC);
   
   // We don't have wisdom, so compute a plan
   if ( NULL == p && opt->plan_rigor > 1 ) {
      // Make copy of A for planning (this doesn't actually occur for FFTW_ESTIMATE, so no copy)
      if (DEBUG_OR_VERB(1))
         DEBUG_PRINTF("Making copy of A to create FFTW plan.\n",NULL);
      clock_t end,begin=clock();
      double *A_cpy = (double *)mxMalloc(m*n*sizeof(*A));

      // Compute the plan
      if (DEBUG_OR_VERB(1))
         DEBUG_PRINTF("Starting FFTW planning.\n",NULL);
      p = fftw_plan_many_r2r(1, trans_lens, (int)n, A_cpy, NULL, 1, (int)m, A_cpy, NULL, 1, (int)m,
		   kind, fftw_flags);
      if (DEBUG_OR_VERB(2))
         DEBUG_PRINTF("fftw_alignment_of(A) = %d; fftw_alignment_of(A_cpy) = %d.\n",
               fftw_alignment_of(A), fftw_alignment_of(A_cpy));
      end = clock();
      if (DEBUG_OR_VERB(1))
         DEBUG_PRINTF("FFTW planning time = %e.\n", (double)(end-begin)/CLOCKS_PER_SEC);
      mxFree(A_cpy);
   }
   else if ( NULL == p && opt->plan_rigor == 1 ) {
      // Compute the plan (using FFTW_ESTIMATE so no need to copy)
      clock_t end,begin=clock();
      p = fftw_plan_many_r2r(1, trans_lens, (int)n, A, NULL, 1, (int)m, A, NULL, 1, (int)m,
		   kind, fftw_flags);
      end = clock();
      if (DEBUG_OR_VERB(1))
         DEBUG_PRINTF("FFTW planning time = %e.\n", (double)(end-begin)/CLOCKS_PER_SEC);
   }
   else
      if (DEBUG_OR_VERB(1))
         DEBUG_PRINTF("Suitable wisdom found to create plan\n",NULL);


   if ( NULL == p) {
      mex_error("error making FFTW plan");
   }
  
   if ( opt->wisdom_only ) { 
      if (DEBUG_OR_VERB(1))
         DEBUG_PRINTF("Only computing wisdom; returning.\n",NULL);
      return;
   }

   // Execute the plan
	for ( it=0; it<n_its; ++it ) {
      
      // Rescale to get desired normalization (orthogonal transform; matches MATLAB)
      if ( nrm_before ) {
         #if __USE_OMP__
         #pragma omp parallel for num_threads((opt->nthreads)) private(i,mj)
         #endif
         for (j=0; j<n; ++j) {
            mj=m*j;
            A[mj] *= nrm_1;
            for (i=1; i<m-1; ++i)
               A[mj+i] *= nrm_mid;
            A[mj+m-1] *= nrm_end;
         }
      }
      
      // Fourier-like goes before D
      if ( !darr_before_F ) { 
		   // Apply transform to columns
		   fftw_execute_r2r(p, A, A);
      }
 
      if ( opt->use_diags ) {
		   // Pointer to proper column of darr array
         if ( !darr_reversed ) // forward through darr
            darr_it = darr + it*m;
         else
            darr_it = darr + (n_its-1-it)*m; // backward through darr
		   
         // Apply +-1 to rows
         #if __USE_OMP__
         #pragma omp parallel for num_threads((opt->nthreads)) private(i,mj)
         #endif
		   for ( j=0; j<n; ++j ) {
            mj = m*j;
            for ( i=0; i<m; ++i ) {
               A[mj+i] *= darr_it[i];
            }
		   }
      }
      
      // Fourier-like goes after D
      if ( darr_before_F ) {
		   // Apply transform to columns
		   fftw_execute_r2r(p, A, A);
      }
   
      // Rescale to get desired normalization (orthogonal transform; matches MATLAB)
      if ( !nrm_before ) {
         #if __USE_OMP__
         #pragma omp parallel for num_threads((opt->nthreads)) private(i,mj)
         #endif
         for (j=0; j<n; ++j) {
            mj=m*j;
            A[mj] *= nrm_1;
            for (i=1; i<m-1; ++i)
               A[mj+i] *= nrm_mid;
            A[mj+m-1] *= nrm_end;
         }
      }
   }

	fftw_destroy_plan(p);
}

void mix_cols_blocked(double *A, double *darr, size_t m, size_t n, size_t n_its, mix_opt *opt) {
		
   size_t bs = opt->block_size;
   size_t n_blocks, block_rem;

   if ( bs > n ) {
      if (DEBUG_OR_VERB(1))
         DEBUG_PRINTF("block size, %d, greater than number of columns, %d.  "
               "Continuing with block_size=%d.\n", bs, n, n);
      opt->block_size=n; // we call the monotlithic code, so this doesn't do anything
      mix_cols(A, darr, m, n, n_its, opt);
      return;
   }

   n_blocks = n/bs;
   block_rem = n%bs;

   // Import wisdom file (if provided by user)
   import_wisdom(opt);
   
   //TODO: investigate better ways to use OMP
   for (size_t b=0; b < n_blocks; ++b)
      mix_cols(A+b*bs*m, darr, m, bs, n_its, opt);
   if (block_rem)
      mix_cols(A+n_blocks*bs*m, darr, m, block_rem, n_its, opt);
   
   // Save wisdom to file (if requested by user)
   export_wisdom(opt);
}


//TODO: This routine works in place on A; when is it better to transpose back+forth?
void mix_rows_inplace(double *A, double *darr, size_t m, size_t n, size_t n_its, mix_opt *opt) {
	mex_error("mix_rows is not implemented.");
}


mxArray* default_options_struct(mix_opt *opt) {
   
   const char *fieldnames[] = {"transpose", "side", "transform",
      "plan_rigor", "wisdom_file", "wisdom_only",
      "use_diags", "verbosity", "nthreads", "block_size"};
   
   mxArray *optmx = mxCreateStructMatrix(1,1,10,fieldnames);
   
   mxSetField(optmx,0, "transpose", mxCreateLogicalScalar(opt->transpose));
   mxSetField(optmx,0, "side", mxCreateString(opt->side));
   mxSetField(optmx,0, "transform", mxCreateString(opt->transform));

   mxSetField(optmx,0, "plan_rigor", mxCreateDoubleScalar(opt->plan_rigor));
   mxSetField(optmx,0, "wisdom_file", mxCreateString(opt->wisdom_file));
   mxSetField(optmx,0, "wisdom_only", mxCreateLogicalScalar(opt->wisdom_only));

   mxSetField(optmx,0, "use_diags", mxCreateLogicalScalar(opt->use_diags));
   mxSetField(optmx,0, "verbosity", mxCreateDoubleScalar(opt->verbosity));
   mxSetField(optmx,0, "nthreads", mxCreateDoubleScalar(opt->nthreads));
   mxSetField(optmx,0, "block_size", mxCreateDoubleScalar(opt->block_size));

   return optmx;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   double *darr, *A, *B;
   mxArray *darr_mx;
   size_t m,n,n_its;
   bool generated_darr=false;

   // Default options
   mix_opt opt = {.transpose=false, .side="left", .transform="dct",
      .plan_rigor=1, .wisdom_file="", .wisdom_only=false,
		.use_diags=true, .verbosity=0, .nthreads=1, .block_size=250};

   /****************************/
   /* Check inputs and outputs */
   /****************************/
   // Check A
   if ( nrhs == 0 ) { // return default options struct
      plhs[0] = default_options_struct(&opt);
      return;
   }
   mxClassID class = mxGetClassID(prhs[0]);
   if ( !(class == mxDOUBLE_CLASS || class == mxSINGLE_CLASS) )
      mex_error("we only handle real single and double arrays");
   if ( mxIsEmpty(prhs[0]) ) mex_error("A is empty");
   if ( mxGetNumberOfDimensions(prhs[0]) != 2 ) mex_error("A should be an m x n matrix");
   
   A = mxGetPr(prhs[0]);
   m = mxGetM(prhs[0]);
   n = mxGetN(prhs[0]);

   // Read in options struct
	// {{{
   if ( nrhs >= 3 ) { 
      mxArray *tmp;

      // transpose
      tmp = mxGetField(prhs[2], 0, "transpose");
      if ( tmp ) {
         if ( !mxIsScalar(tmp) )
            mex_error("opt.transpose should be either 0 or 1");

         if ( mxGetScalar(tmp) == 0. )
            opt.transpose = false;
         else if ( mxGetScalar(tmp) == 1. )
            opt.transpose = true;
         else
            mex_error("opt.transpose should be either 0 or 1");
      }

      // side
      tmp = mxGetField(prhs[2], 0, "side");
      if ( tmp ) {
         opt.side = mxArrayToString(tmp);
         if ( !(opt.side) ) mexErrMsgTxt("error converting opt.side to string.");
         for (char *p=opt.side;*p;++p) *p=tolower(*p);
         if ( !(strncmp(opt.side,"l",2)==0 || strncmp(opt.side,"left",5)==0 ||
                strncmp(opt.side,"r",2)==0 || strncmp(opt.side,"right",6)==0) )
            mex_error("opt.side should be 'l'/'left' or 'r'/'right'");
      }

      // transform
      tmp = mxGetField(prhs[2], 0, "transform");
      if ( tmp ) {
         opt.transform = mxArrayToString(tmp);
         if ( !(opt.transform) ) mexErrMsgTxt("error converting opt.transform to string.");
         for (char *p=opt.transform;*p;++p) *p=tolower(*p);
         if ( !(strncmp(opt.transform,"dct",4)==0 || strncmp(opt.transform,"idct",5)==0 ) )
            mex_error("opt.transform should be 'DCT' or 'IDCT'");
      }

      // plan_rigor
      tmp = mxGetField(prhs[2], 0, "plan_rigor");
      if ( tmp ) {
         if ( !mxIsScalar(tmp) )
            mex_error("opt.plan_rigor should be 1, 2, 3, or 4");
         
         int pr = (int)mxGetScalar(tmp);
         if ( pr >= 1 && pr <= 4 ) 
            opt.plan_rigor = pr;
         else
            mex_error("opt.plan_rigor should be 1, 2, 3, or 4");
      }

      // wisdom_file
      tmp = mxGetField(prhs[2], 0, "wisdom_file");
      if ( tmp ) {
         if ( !mxIsChar(tmp) )
            mexErrMsgTxt("error converting opt.wisdom_file to string.");
         opt.wisdom_file = mxArrayToString(tmp);
         if ( !(opt.wisdom_file) ) mexErrMsgTxt("error converting opt.wisdom_file to string.");
      }

      // wisdom_only
      tmp = mxGetField(prhs[2], 0, "wisdom_only");
      if ( tmp ) {
         if ( !mxIsScalar(tmp) )
            mex_error("opt.wisdom_only should be either 0 or 1");

         if ( mxGetScalar(tmp) == 0. )
            opt.wisdom_only = false;
         else if ( mxGetScalar(tmp) == 1. )
            opt.wisdom_only = true;
         else
            mex_error("opt.wisdom_only should be either 0 or 1");
      }
      
      // use_diags
      tmp = mxGetField(prhs[2], 0, "use_diags");
      if ( tmp ) {
         if ( !mxIsScalar(tmp) )
            mex_error("opt.use_diags should be either 0 or 1");

         if ( mxGetScalar(tmp) == 0. )
            opt.use_diags = false;
         else if ( mxGetScalar(tmp) == 1. )
            opt.use_diags = true;
         else
            mex_error("opt.use_diags should be either 0 or 1");
      }

      // verbosity
      tmp = mxGetField(prhs[2], 0, "verbosity");
      if ( tmp ) {
         if ( !mxIsScalar(tmp) )
            mex_error("opt.verbosity should be 0, 1, or 2");
         
         int vrb = (int)mxGetScalar(tmp);
         if ( vrb >= 0 && vrb <= 2 ) 
            opt.verbosity = vrb;
         else
            mex_error("opt.verbosity should be 0, 1, or 2");
      }
      
      // nthreads
      tmp = mxGetField(prhs[2], 0, "nthreads");
      if ( tmp ) {
         if ( !mxIsScalar(tmp) )
            mex_error("opt.nthreads should be an integer.");
         
         int nt = (int)mxGetScalar(tmp);
         if ( nt >= -1) 
            opt.nthreads = nt;
         else
            mex_error("opt.nthreads should be an integer.");
      }

      // block_size
      tmp = mxGetField(prhs[2], 0, "block_size");
      if ( tmp ) {
         if ( !mxIsScalar(tmp) )
            mex_error("opt.block_size should be an integer.");
         
         size_t bs = (size_t)mxGetScalar(tmp);
         if ( bs >= 0) 
            opt.block_size = bs;
         else
            mex_error("opt.block_size should be an integer.");
      }

	}

   if (DEBUG >= 1 || opt.verbosity >= 1)
      DEBUG_PRINTF("\ninterface options:\n -> transpose=%d\n -> side=\"%s\"\n"
            " -> transform=\"%s\"\n -> plan_rigor=%d\n"
            " -> wisdom_file=\"%s\"\n -> wisdom_only=%d\n\n -> use_diags=%d\n"
            " -> verbosity=%d\n -> nthreads=%d\n -> block_size=%d\n",
            opt.transpose, opt.side, opt.transform, 
            opt.plan_rigor, opt.wisdom_file, opt.wisdom_only,
            opt.use_diags, opt.verbosity, opt.nthreads, opt.block_size);
	// }}}
	
   // n_its+generate darr or get pointer to darr's data
   bool is_scalar = false;
   if ( nrhs >= 2 ) {
      if ( mxIsEmpty(prhs[1]) ) // if darr=[], don't apply darr
         opt.use_diags = false;
      is_scalar = mxIsScalar(prhs[1]);
   }
   if ( nrhs < 2 || is_scalar ) { // not specified or n_its given
      if ( is_scalar ) { // n_its is provided
         n_its = (size_t)mxGetScalar(prhs[1]);
         if ( (double)n_its != mxGetScalar(prhs[1]) )
            mex_error("n_its should be an integer.");
      }
      else // default value of n_its
         n_its = 1;
      
      // Use MATLAB's `rand()` to populate darr with Uniform[0,1] 
      mxArray* mcm_prhs[2];
      size_t trans_dim;
      if ( *opt.side == 'l' )
         trans_dim = m;
      else
         trans_dim = n;
      
      mcm_prhs[0] = mxCreateDoubleScalar(trans_dim);
      mcm_prhs[1] = mxCreateDoubleScalar(n_its);
      mexCallMATLAB(1, &darr_mx, 2, mcm_prhs, "rand");
      darr = mxGetPr(darr_mx);
      
      // Threshold to get +- 1 values
      for (unsigned i=0; i < trans_dim*n_its; ++i) {
         darr[i] = (darr[i] >= 0.5) ? 1. : -1.;
      }

      generated_darr = true;
   }
   else { // Use the given darr
      darr = mxGetPr(prhs[1]);
      n_its = mxGetN(prhs[1]);
   
      if ( *opt.side == 'r' && n != mxGetM(prhs[1]) )
         mex_error("dimension 2 of A and dimension 1 of darr should be equal.");
      else if ( *opt.side == 'l' && m != mxGetM(prhs[1]) )
         mex_error("dimension 1 of A and dimension 1 of darr should be equal.");
   }

   // Set output pointers
   if ( nlhs >= 2 ) {
      if ( !generated_darr )
         mex_error("darr output specified but is only returned when generating "
               "a new darr.");
      plhs[1] = darr_mx;
   }

   /***************/
   /* Do the work */
   /***************/
   // Set up FFTW threads
   int s = fftw_init_threads();
   if ( s == 0 ) 
      mex_error("Error setting up FFTW threads\n.");
   fftw_plan_with_nthreads(opt.nthreads);
   if (DEBUG >= 1 || opt.verbosity >= 1)
      DEBUG_PRINTF("FFTW planning with %d threads\n", opt.nthreads);
   
   // Call the worker routines
   if ( *opt.side == 'l' ) {
      // Make a copy of A to work on
      plhs[0] = mxDuplicateArray(prhs[0]); B = mxGetPr(plhs[0]);
      
      // mix the columns
      if ( opt.block_size > 0 )
         mix_cols_blocked(B, darr, m, n, n_its, &opt);
      else
         mix_cols(B, darr, m, n, n_its, &opt);
   }
   else {
      // Make a copy of the transpose of A to work on
      mxArray *At_mx;
      mexCallMATLAB(1, &At_mx, 1, (mxArray **) &prhs[0], "transpose");
      double *At = mxGetPr(At_mx);

      // Mix the columns 
      if ( opt.block_size > 0 )
         mix_cols_blocked(At, darr, n, m, n_its, &opt);
      else
         mix_cols(At, darr, n, m, n_its, &opt);

      // Initialize output array and copy over transpose of mixed At
      mexCallMATLAB(1, &plhs[0], 1, &At_mx, "transpose");

   }
   
   /************/
   /* Clean up */
   /************/
   fftw_cleanup_threads();
   //TODO: only free if we called mxArrayToString
   //mxFree(opt.side);
   //mxFree(opt.transform);
   //mxFree(opt.wisdom_file);

   return;
}

