function A = extkahan(l, phi, mu)
% EXTKAHAN  Extended Kahan matrix
%
% Inputs: l     size of three block matrices; returned matrix size is [3l 3l]
%               l, l/12, or l/20 should be a power of two.
%         phi   angle parameter [zeta^2+phi^2 == 1 in [3] and
%               phi > (4l-1)^(-1/4) ]
%               defaults to 0.285
%         mu    perturbation parameter (0 < mu << 1)
%               defaults to 20*eps(1)/sqrt(3l)
%
%  [3] - Efficient Algorithms for Computing a Strong Rank-Revealing QR Factorization,
%        M. Gu and S. Eisenstat, SIAM J. on Sci. Comp., 1996

if nargin < 1, error('Not enough input arguments.'); end
%if l ~= bitshift(1, round(log2(l))), error('l should be a power of 2.'); end
[f,e] = log2([l l/12 l/20]);
k = find(f == 1/2 & e > 0);
if isempty(k), error('l, l/12, or l/20 should be a power of two.'); end
if nargin < 2, phi = 0.285; end
if nargin < 3, mu = 20*eps(1)/sqrt(3*l); end

if ( phi > 1 || phi < (4*l-1)^(-.25) )
   warning('Paramter phi should be (4l-1)^(-.25) < phi <= 1.');
end

zeta = sqrt(1-phi^2);

H = hadamard(l);

A = diag(zeta.^(0:3*l-1)) * ...
    [eye(l) -phi*H zeros(l); zeros(l) eye(l) phi*H; zeros([l 2*l]) mu*eye(l)];

end
