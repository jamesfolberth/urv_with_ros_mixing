function [A,s] = devil(n, l, alpha)
% DEVIL  Devil's staircase matrix of size n with step length l
%
% Inputs: n      size of matrix
%         l      stair length (default: 20)
%         alpha  jump factor 10^-alpha (default: 0.6)
%
% Outputs: A   matrix
%          s   vector of true singular values
%
%  [2] - Communication Avoiding Rank Revealing QR with column pivoting,
%        J. Demmel, L. Grigori, M. Gu, et al., SIAM J. Matrix Anal. Appl., 2015

if nargin < 1, error('Not enough input arguments.'); end
if nargin < 2, l = 20; end
if nargin < 3, alpha = 0.6; end

s = zeros(n,1);
n_stairs = floor(n/l);
for i=1:n_stairs
   s(1+l*(i-1):l*i) = -alpha*(i-1);
end
if n_stairs > 0
   s(l*n_stairs:end) = -alpha*(n_stairs);
end
s = 10.^s;
A = orth(rand(n)) * diag(s) * orth(randn(n));

end
