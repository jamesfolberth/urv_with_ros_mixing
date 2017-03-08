function [A,s] = stewart2(N,smin)
% STEWART2  Random m x n matrix of rank min(m,n)//2 with decaying svals and noise from [7]. 
%
% Return the matrix A = U*S*V + 0.1*smin*rand(m,n), where U,V are Haar random 
% orthogonal matrices, and s is constructed as follows:
%   p = min(m,n); nhalf = floor(p/2);
%   s = zeros([p 1]); s(1:nhalf) = fliplr(logspace(log10(smin), log10(1), nhalf));
%
% Inputs: N     size of matrix (scalar or [m n])
%         smin  smallest non-zero singular value (default: 1e-3)
%
% Outputs: A   matrix
%          s   vector of true singular values
%

if nargin < 2, smin = 1e-3; end

m = N(1); n = N(length(N));
p = min(m,n);
[U,~] = qr(randn(m,p),0);
[V,~] = qr(randn(n,p),0);
nhalf = floor(p/2);
s = zeros([p 1]); s(1:nhalf) = fliplr(logspace(log10(smin), log10(1), nhalf));
A = U*diag(s)*V' + 0.1*smin*rand(m,n);

end
