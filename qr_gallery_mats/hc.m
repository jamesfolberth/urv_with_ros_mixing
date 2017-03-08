function [A,s] = hc(N)
% HC  Random m x n matrix with prescribed singular values
%
% Return the matrix A = U*S*V, where U,V are Haar random orthogonal matrices
% and s(1) = 100; s(2) = 10; s(3:end) = logspace(-2,-8,p-2).
%
% Inputs: N     size of matrix (scalar or [m n])
%
% Outputs: A   matrix
%          s   vector of true singular values
%

m = N(1); n = N(length(N));
p = min(m,n);
[U,~] = qr(randn(m,p),0);
[V,~] = qr(randn(n,p),0);
s = zeros([n 1]);
if p >= 1, s(1) = 100; end
if p >= 2, s(2) = 10; end
if p >= 3, s(3:end) = logspace(-2,-8,p-2); end
A = U*diag(s)*V';

end

