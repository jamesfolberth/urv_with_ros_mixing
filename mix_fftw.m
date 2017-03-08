function [varargout] = mix_fftw(A, D_diags, opt)
% MIX_FFTW  Performs mixing with transforms from FFTW and Rachemacher +- 1
%
% MIX_FFTW computes the operation A*V, A*V', V*A, or V'*A where
%
%   V = (F*D_k)*...*(F*D_1),
%
% F is a fast orthogonal transform from FFTW specified by opt.transform,
% and D_i is a diagonal matrix of random +- 1.  The side on which V is applied 
% is specified by opt.side; the use of V or V' is specified by opt.transpose;
% the transform used in V is specified by opt.transform.  The diagonals of the
% D_i are stored in D_diags.
%
% Inputs: A                   m x n real matrix
%         D_diags or n_its    m x n_its matrix of random +-1 (n x n_its for opt.side='right') or
%                             number of iterations (when generating D_diags) or
%                             pass in empty matrix [] to skip using +-1 in mixing
%
%         opt              options structure
%           opt.transpose     1 to use V'; 0 to use V
%           opt.side          'left' of A or 'right' of A
%           opt.transform     'DCT', 'IDCT' (transform used by V)
%
%           opt.plan_rigor    1,2,3,4.  Higher means better runtime but (much) longer planning.
%           opt.wisdom_file   filepath to file with FFTW wisdom (will create if it doesn't exist)
%           opt.wisdom_only   If 1, compute plan, cache wisdom in opt.wisdom_file, and return 
%
%           opt.use_diags     skip using D_diags in the mixing (will still do n_its iterations)
%           opt.verbosity     0,1,2.  Higher means more verbose output.  0 should be quiet.
%           opt.nthreads      number of threads to use in FFTW and OpenMP.
%           opt.block_size    perform the transform on column blocks of this size
%                             set opt.block_size=0 to use the full matrix.
%
%     Both D_diags and n_its are optional:
%        mix_rows_dct(A) sets n_its=1 and generates D_diags
%        mix_row_dct(A,n_its) generates D_diags with n_its
%        mix_rows_dct(A,D_diags) uses the given D_diags with n_its=size(D_diags,2)
%
%     If D_diags is generated (i.e. not input), it is returned in the 2nd output argument
%
%     If no arguments are given, a struct with the default options is returned
%
% Outputs: B     mixed A
%          D_diags  output only if generated
%
%     Note that if opt.wisdom_only == 1, neither B nor D_diags are output
%
%     If no input arguments are given, the default options structure is returned, i.e.
%     >> mix_opt = mix_fftw();
%
% Examples:
%
%   %% Use randomized URV with ROS mixing to solve LS problems
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Generate a random data matrix and RHS vectors:
%   >> A = randn(500,50);
%   >> b = randn(500,3);
%
%   % Perform ROS mixing:
%   >> mix_opt = mix_fftw();                       % get the default options
%   >> mix_opt.side = 'right';                     % want to mix as A_hat = A*V';
%   >> mix_opt.transpose = 1;
%   >> mix_opt.transform = 'DCT';                  % the transform in V is a DCT
%   >> [A_mix, D_diags] = mix_fftw(A, 2, mix_opt); % 2 ROS iterations
%
%   % Unpivoted QR factorization (skipping pre-sort) and solve with backward substitution:
%   >> [U,R] = qr(A_mix,0);
%   >> y_mix = backsub(R, U.'*b); % backsub.m performs backward subs. in MATLAB
%
%   % Unmix `y_mix` to find the solution `x=V^T*y_mix`:
%   >> x = mix_fftw(y_mix, D_diags, unmix_opt); % x = V'*y
%
%   % Check the residual:
%   >> norm(A*x-b, 'fro')                                 
%   ans =
%      36.2982
%   >> norm(A*(A\b)-b, 'fro') % should produce the same residual
%   ans =
%      36.2982
%
%
%   %% Cache wisdom
%   %%%%%%%%%%%%%%%
%   >> opt_mix.plan_rigor = 4;                        
%   >> opt_mix.wisdom_file = 'wisdom.fftw';           % will create the file, but not dirs
%   >> [A_mix, D_diags] = mix_fftw(A, 2, opt_mix);    % FFTW planning takes a bit; caches wisdom
%   >> [A_mix, D_diags] = mix_fftw(A, 2, opt_mix);    % loads wisdom file; planning is very fast
%
