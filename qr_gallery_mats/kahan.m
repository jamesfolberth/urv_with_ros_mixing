function A = kahan(N, theta, tau)
% KAHAN  Kahan matrix with Demmel/Grigori perturbation
%
% Inputs: N     size of matrix (scalar or [m n])
%         theta angle parameter (c^2+s^2 == 1 in [1])
%               defaults to acos(0.6)
%         tau   perturbation parameter (different than MATLAB's 'kahan')
%               defaults to 1e-7, for which QRCP doesn't pivot for at least N=100
%               change to 1e-15, and QRCP should pivot.
%
%  [1] - Communication Avoiding Rank Revealing QR, J. Demmel, L. Grigori,
%        M. Gu, et al., SIAM J. Matrix Anal. Appl., 2015

if nargin < 1, error('Not enough input arguments.'); end
if nargin < 2, theta = acos(0.6); end
if nargin < 3, tau = 1e-7; end

m = N(1);
n = N(length(N));

c = cos(theta); s = sin(theta);

A = eye(n) - c*triu(ones(n), 1);
%A = diag(s.^(0:n-1))*A*diag((1-tau).^(0:n-1));
A = bsxfun(@times, bsxfun(@times,A,s.^(0:n-1).'), (1-tau).^(0:n-1)); % probably less roundoff
%norm(A - bsxfun(@times, bsxfun(@times,A_,s.^(0:n-1).'), (1-tau).^(0:n-1)), 'fro')

end
