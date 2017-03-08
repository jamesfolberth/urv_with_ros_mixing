function [A,s] = break1(N)
% BREAK1  Random m x n matrix with prescribed singular values
%
% Return the matrix A = U*S*V, where U,V are Haar random orthogonal matrices
% and S = diag([ones([p-1,1]); 1e-9]), where p=min(m,n).
%
% Inputs: N     size of matrix (scalar or [m n])
%
% Outputs: A   matrix
%          s   vector of true singular values
%

m = N(1); n = N(length(N));
p = min(m,n);
if p < 2, error('min(m,n) should be >= 2.'); end
[U,~] = qr(randn(m,p),0);
[V,~] = qr(randn(n,p),0);
s = ones([p,1]);
s(end) = 1e-9;
A = U*diag(s)*V';

end

