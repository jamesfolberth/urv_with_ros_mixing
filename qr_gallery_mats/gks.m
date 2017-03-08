function A = gks(n)
% GKS  Upper-triangular matrix with specified values
%
% A is an upper-triangular matrix where
%  1. the jth diagonal element is 1/sqrt(j)
%  2. the (i,j) element is -1/sqrt(j) for j>i (zero otherwise)
%
% Inputs: n     size of n x n matrix

A = diag(1./sqrt(1:n)) - triu(bsxfun(@times, ones(n), 1./sqrt(1:n)),1);

end
