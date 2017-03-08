# Randomized URV Factorization with ROS Mixing
Stephen Becker, James Folberth, Laura Grigori

https://arxiv.org/abs/1703.02499

## Implementation of ROS Mixing
We have implemented random orthogonal system (ROS) mixing in the function `mix_fftw`.  `mix_fftw` is written in C with a MATLAB MEX interface.  We use [FFTW3](http://www.fftw.org/) for fast orthogonal transforms (e.g., discrete cosine/sine transform, discrete Hartley transform), and (GNU) [OpenMP](http://www.openmp.org/) for additional parallelization.  For improved speed, `mix_fftw` optionally caches FFTW wisdom to a file; see `help mix_fftw` for more information and usage.

We have a simple script, `setup_mex`, which will attempt to build our MEX files and download dependencies for our other experiments.  The interface to `mix_fftw` and its various options are described in `mix_fftw.m`; the help text can be viewed in MATLAB with `help mix_fftw`.

NOTE: Currently, we have only implemented mixing with the DCT.  If there is interest, we are happy to implement other transforms.

## Example Usage - Least Squares
As an illustrative example, we show how to use `RURV_ROS` to solve least-squares problems.  We use `ls_ROS.m` in our testing, which uses `mix_fftw` for mixing and `qr_ls`, a MEX function to solve least-squares with unpivoted QR.

Consider solving the overdetermined least-squares problem `min norm(A*x-b,'fro')` with multiple right-hand side vectors in `b`.  `RURV_ROS` produces a factorization of the form `A=U*R*V`.  We compute the mixed matrix `A_mix = A*V^T` and (implicitly) make the substitution `y_mix=V*x`.  We then solve for `y_mix` with `y_mix = R^(-1)*U^T*b`, where `A_mix = U*R` is computed with an unpivoted QR and `R^(-1)` is applied with backward substitution.

Generate a random matrix and three right-hand side vectors:
```matlab
A = randn(500,50);
b = randn(500,3);
```
Perform ROS mixing with `mix_fftw`:
```matlab
mix_opt = mix_fftw();                       % get the default options
mix_opt.side = 'right';                     % want to mix as A_hat = A*V';
mix_opt.transpose = 1;
mix_opt.transform = 'DCT';                  % the transform in V is a DCT
[A_mix, D_diags] = mix_fftw(A, 2, mix_opt); % 2 ROS iterations
```
Unpivoted QR factorization (skipping pre-sort) and solve with backward substitution:
```matlab
[U,R] = qr(A_mix,0);
y_mix = backsub(R, U.'*b); % backsub.m performs backward subs. in MATLAB
```
Unmix `y_mix` to find the solution `x=V^T*y_mix`:
```matlab
x = mix_fftw(y_mix, D_diags, unmix_opt); % x = V'*y
```
Check the residual:
```
>> norm(A*x-b, 'fro')                                 
ans =
   36.2982
>> norm(A*(A\b)-b, 'fro') % should produce the same residual
ans =
   36.2982
```

## Our Drivers
* `ls_ROS.m`

   Driver to solve least-squares problems with `RURV_ROS` (or `RVLU_ROS` if seeking the minimum norm solution to an underdetermined system).  This function uses `mix_fftw` and `qr_ls` for ROS mixing and solving LS problems with unpivoted QR, respectively.
* `ls_test.m`

   Driver for various least-squares tests.  We experiment with runtimes and computing accurate solutions to underdetermined systems with correlated columns.
* `qr_timing.m`

   Driver for timing experiments with LAPACK's unpivoted QR, column pivoted QR, and `RURV_ROS`.  We use `qr_fact_only` to call LAPACK's `dgeqrf` or `dgeqp3` without building the orthonormal factor, which costs the same for all three factorizations.
* `rr_test.m`

   Driver to experiment with the rank-revealing conditions and the behavior of `RURV_ROS` compared to `RURV_Haar`.
* `rval_ratio_test.m`

   How well do the R-values (diagonal elements of the `R` factor) approximate the SVD.  The answer is "not well" for `RURV_ROS`, but this example illustrates the effect of using two mixing steps instead of the usual one mixing step.
* `qlp_test.m`

   Experiments with Stewart's QLP approximation to the SVD.
* `column_norm_smoothing.m`

   Visualize the effect of mixing on the column norms.

## `qr_gallery`
In our research, we saw a variety of matrices that others have tested with.  We collected a number of these matrices and provide an interface to them called `qr_gallery`.  The syntax is similar to MATLAB's built-in `gallery` function.  See the help text with `help qr_gallery` in MATLAB.

Some of the matrices accessible from `qr_gallery` are from P.C. Hansen's [Regularization Tools](http://www.imm.dtu.dk/~pcha/Regutools/).  They are included here in `qr_gallery_mats/` with the appropriate attribution and license.

## External Code and Dependencies.
In our least-squares tests, we compare against [BLENDENPIK](http://epubs.siam.org/doi/abs/10.1137/090767911).  MATLAB MEX code is available [here](https://www.mathworks.com/matlabcentral/fileexchange/25241-blendenpik).

We print our figures in MATLAB with [export_fig](https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig).

We use `DGEQPX`, a rank-revealing QR, with a MATLAB MEX interface from Foster and Liu's rank-revealing factorization comparison paper.  Their paper and code are available [here](http://www.math.sjsu.edu/~foster/rankrevealingcode.html).  We attempt to automatically download and install their code in the function `add_rrqr_comparison`, which is called by `setup_mex`.

We use a MEX interface to MATLAB's LAPACK from [here](https://www.mathworks.com/matlabcentral/fileexchange/16777-lapack).  We attempt to automatically download and install the interface in the function `add_lapack`, which is called by `setup_mex`.

We occasionally use [`suptitle.m`](https://www.mathworks.com/matlabcentral/fileexchange/45049-bland-altman-and-correlation-plot/content/BlandAltman/suptitle.m) to make the title when using a figure with subplots.
