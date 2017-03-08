function [A,s] = condk(N, kappa)
% CONDK  Matrix with condition number kappa
%
% Inputs: N     size of matrix (scalar or [m n])
%         kappa condition number
%
% Outputs: A
%          s    singular values

if nargin < 1, error('Not enough input arguments.'); end
if nargin < 2, kappa=1e6; end 

if kappa < 1, error('Condition number should be >= 1.'); end

m = N(1);
n = N(length(N));
p = min(m,n);

%[U,~,V] = svd(randn(m,n),'econ');
% a bit faster ...
[U,~] = qr(randn(m,p),0);
[V,~] = qr(randn(n,p),0);
s = linspace(1/sqrt(kappa), sqrt(kappa), p); s = s(end:-1:1);
S = diag(s);
A = U(:,1:p)*S*V(:,1:p)';

end
