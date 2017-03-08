function A = stewart(mu)
% STEWART  Small example matrix from G.W. Stewart
%
% Inputs: mu Parameter for two large rows.
%            defaults to 1e12
%
%  [4] - Matrix Algorithms. Volume I: Basic Decompositions.
%        G.W. Stewart, SIAM, 1998

if nargin < 1, mu = 1e12; end

A = [1 1 1;
     1 3 1;
     1 -1 1;
     1 1 1;
     mu mu mu;
     mu mu -mu];

end
