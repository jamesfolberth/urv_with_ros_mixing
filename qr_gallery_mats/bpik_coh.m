function [varargout] = bpik_coh(N, kappa, coh, pert)
% BPIK_COH  In-/Semi-/-coherent matrices used to test BLENDENPIK
%
% Inputs: N     size of matrix (scalar or [m n])
%         kappa approximate condition number
%         coh   'i' - incoherent
%               's' - semicoherent
%               'c' - coherent
%         pert  perturbation of pert*ones(N)
%               (optional: default 1e-8)
%
%   These matrices are usually used for overdetermined systems (i.e. m > n).
%   The coherent mode (coh=='c') will work for m < n.
%
% Outputs: A
%          s    singular values (only if coh=='i')
% 
%  [5] - BLENDENPIK: supercharging LAPACK's least-squares solver.
%        H. Avron, et al., SIAM J. on Sci. Comp., 2010

if nargin < 1, error('Not enough input arguments.'); end
if nargin < 2, kappa = 1e6; end
if nargin < 3, coh = 'i'; end
if nargin < 4, pert = 1e-8; end

if kappa < 1, error('Condition number should be >= 1.'); end
if strcmp(coh,'')==1, error('coh should be either ''i'', ''s'', or ''c'''); end

m = N(1);
n = N(length(N));
p = min(m,n);

if strncmp(coh,'i',1) == 1 % incoherent
   [U,~] = qr(randn(m,p),0);
   [V,~] = qr(randn(n,p),0);
   s = linspace(1/sqrt(kappa), sqrt(kappa), p);
   A = U(:,1:p)*diag(s)*V(:,1:p)';
   varargout{1} = A;
   if nargout > 1, varargout{2} = s; end

elseif strncmp(coh,'s',1) == 1 % semicoherent
   if m < n
      error('We don''t support underdetermined (m<n) systems for the semicoherent case');
   end
   m_ = floor(m - n/2);%TODO we could let the user give the factor 1/2 in n/2
   n_ = floor(n/2);
   p_ = min(m_,n_);
   [U,~] = qr(randn(m_,p_),0);
   [V,~] = qr(randn(n_,p_),0);
   s = linspace(1/sqrt(kappa), sqrt(kappa), p_);
   B = U(:,1:p_)*diag(s)*V(:,1:p_)';
   A = [B zeros(m_, n-n_);zeros(m-m_,n_) eye(m-m_,n-n_)];
   A = A + pert;
   varargout{1} = A;

elseif strncmp(coh,'c',1) == 1 % coherent
   %TODO we don't use the kappa argument
   if m < n
      error('We don''t support underdetermined (m<n) systems for the coherent case');
   end
   D = diag(randn(n,1));
   A = [D;zeros(m-n,n)];
   A = A + pert;
   varargout{1} = A;

else
   error('Unrecognized coherence parameter %s', coh);
end

end
