function [A,s] = rankk(N,k,drop)
% RANKK  Random (real Haar) m x n matrix of approximate rank k.
%
% Inputs: N     size of matrix (scalar or [m n])
%         k    numerical rank.  The singular values will decay slowly for 1:k,
%              then will have steep drop and decay slowly again for k+1:end.
%         drop depth of jump in singular values (default: 1e-10) 
%
% Outputs: A   matrix
%          s   vector of true singular values
%

if nargin < 2, error('Not enough input arguments.'); end
if nargin < 3, drop=1e-10; end

m = N(1); n = N(length(N));
p = min(m,n);
if p < k, warning('Requested rank of A should be <= min(m,n)'); k=p; end
[U,~] = qr(randn(m,p),0);
[V,~] = qr(randn(n,p),0);
s = [1./(1:k) drop*1./((k+1:p)-k)].';
A = U*diag(s)*V';

end
