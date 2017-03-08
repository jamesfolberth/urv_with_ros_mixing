% A simple script to build our mex files
% We've succesfully build with GCC 4.4.7 on RHEL
% and GCC 6.3.1 on Arch Linux

% We use the c99 standard only to use the // comment style
mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -O2" qr_ls.c -lmwblas -lmwlapack
mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -O2" qrcp_ls.c -lmwblas -lmwlapack
mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -O2" qr_fact_only.c -lmwblas -lmwlapack

%TODO for more accurate timing experiments, we should link with Matlab's FFTW3
%     However, we had issues with that in development, so we just linked with
%     the system FFTW.  We could instead link with MKL's FFTW, which is probably
%     what Matlab uses internally.
mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -O2 -fopenmp" mix_fftw.c -lgomp -lfftw3_threads -lfftw3 -lm


add_lapack(); % download LAPACK mex interface
add_rrqr_comparison(); % download RRQR codes

