function [varargout] = qr_gallery(matname, varargin)
% QR_GALLERY  Similar to MATLAB's gallery, but specific to testing QRs, URVs, etc.
%
% The syntax is [varargout] = qr_gallery(matname, varargin).  Not all special
% matrices take/return extra arguments.
%
%  baart       - Discretization of 1st kind Fredholm I.E. from [6].
%  bpik_coh    - A random m x n matrix with condition kappa and somewhat 
%                    controllable coherence. Used in [5] to test BLENDENPIK.
%  break1      - Random m x n matrix with prescribed singular values from [2,7].
%  break9      - Random m x n matrix with prescribed singular values from [2,7].
%  condk       - A random m x n matrix with condition number kappa.
%  corrcol     - A random m x n matrix with correlated columns.
%  deriv2      - Computation of 2nd derivative from [6].
%  devil       - Devil's stairs.  Discussed in [2].
%  exponential - Random m x n matrix with exponentially decaying singular values from [2,7].
%  extkahan    - Extended Kahan matrix from [3].
%  foxgood     - Ill-posed test problem of 1st kind Fredholm I.E. from [6].
%  gks         - An upper-triangular matrix specified values from [2].
%  gravity     - 1D gravity surveying problem from [6].
%  hc          - Random m x n matrix with specified singular values from [2].
%  heat        - Inverse heat equation from [6].
%  kahan       - Kahan's example matrix.  Discussed in [1].
%  kahan_rankk - A rank-k modification of Kahan's example.  Modified from 'kahan' [1].
%  phillips    - Phillips' famous test problem from [6].
%  rankk       - A (real Haar) random m x n matrix of approximate rank k.
%  shaw        - 1D image restoration model from [6].
%  spikes      - Test problem with a "spiky" solution from [6].
%  stewart     - An small example due to G.W. Stewart.  Discussed in [4].
%  stewart2    - Random m x n matrix of rank min(m,n)//2 with decaying svals and noise from [7].
%  ursell      - Integral equation with no square integrable solution from [6].
%  wing        - Test problem with a discontinuous solution from [6].
%
% To get more help on special matrices, run qr_gallery() once (to set path)
% and run `help matname`.
%
% References:
%  [1] - Communication Avoiding Rank Revealing QR, J. Demmel, L. Grigori,
%        M. Gu, et al., SIAM J. Matrix Anal. Appl., 2015
%
%  [2] - Communication Avoiding Rank Revealing QR with column pivoting,
%        J. Demmel, L. Grigori, M. Gu, et al., SIAM J. Matrix Anal. Appl., 2015
%  
%  [3] - Efficient Algorithms for Computing a Strong Rank-Revealing QR Factorization,
%        M. Gu and S. Eisenstat, SIAM J. on Sci. Comp., 1996
%
%  [4] - Matrix Algorithms. Volume I: Basic Decompositions.
%        G.W. Stewart, SIAM, 1998
%
%  [5] - BLENDENPIK: supercharging LAPACK's least-squares solver.
%        H. Avron, et al., SIAM J. on Sci. Comp., 2010
%
%  [6] - Regularization Tools version 4.1 (for MATLAB version 7.3).
%        P.C. Hansen, Numerical Algorithms, 46 (2007), pp. 189-194.
%        Code published under a BSD 3-clause license; please see the
%        license in qr_galler_mats/regu_tools_license.txt
% 
%  [7] - A parallel QR factorization algorithm with controlled local pivoting
%        C.H. Bischof, SIAM J. Sci. Stat. Comp., 12 (1991), pp. 36-57.

persistent loaded_path; % init to []
if isempty(loaded_path)
   dir = fileparts(mfilename('fullpath'));
   %fprintf(1, 'Adding %s to path.\n', [dir filesep 'qr_gallery_mats']);
   addpath([dir filesep 'qr_gallery_mats']);
   loaded_path = 1;
end

if nargin == 0, return; end

switch matname
case {'baart', 'bpik_coh', 'break1', 'break9',...
      'condk', 'corrcol',...
      'deriv2', 'devil',...
      'extkahan', 'exponential',...
      'foxgood',...
      'gks', 'gravity',...
      'hc', 'heat',...
      'kahan', 'kahan_rankk',...
      'phillips',...
      'rankk',...
      'shaw', 'spikes', 'stewart', 'stewart2',...
      'ursell',...
      'wing'}
	F = str2func(matname);
	[varargout{1:max(nargout,1)}] = F(varargin{:});
	return

otherwise
	error('Unknown matname: %s', matname)

end

end

