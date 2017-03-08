function [A,s] = exponential(N,alpha)
% EXPONENTIAL  Random m x n matrix with exponentially decaying singular values
%
% Return the matrix A = U*S*V, where U,V are Haar random orthogonal matrices
% and S = diag(alpha.^(0:p-1)), where p=min(m,n).
%
% Inputs: N     size of matrix (scalar or [m n])
%         alpha exponential decay factor (default: 10^(-1/11))
%
% Outputs: A   matrix
%          s   vector of true singular values
%

if nargin < 2, alpha=10^(-1/11); end

m = N(1); n = N(length(N));
p = min(m,n);
[U,~] = qr(randn(m,p),0);
[V,~] = qr(randn(n,p),0);
s = alpha.^(0:p-1).';
A = U*diag(s)*V';

end
