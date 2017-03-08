function [x] = qr_ls(A,b,mode)
% QR_LS  Use unpivoted QR to solve full-rank LS problems
%
% If m >= n, we assume the LS problem is overdetermined.
% Use x = qr_ls(A,b) to solve the LS problem.
%
% If m < n, we assume the LS problem is underdetermined.
% Use x = qr_ls(A,b,mode) to solve the LS problem.
% If mode is 'basic' (default), find the solution of the form [y; zeros].
% If mode is 'minnorm', find the solution x with minimum 2-norm.
%
% Example MATLAB implementation is below.

%   if nargin < 3, mode='basic'; end
%   
%   [m,n] = size(A);
%   if m >= n % overdetermined
%     [Q,R] = qr(A,0); % unpivoted
%     x = backsub(R(1:n,1:n), Q'*b);
%   else % underdetermined
%      if strcmp(mode, 'basic')
%         [Q,R] = qr(A(1:m,1:m),0); % unpivoted
%         x = [backsub(R,Q'*b); zeros(n-m,size(b,2))];
%      elseif strcmp(mode, 'minnorm')
%         [Q,R] = qr(A',0); % unpivoted "LQ" of A
%         x = Q*forwsub(R(1:m,1:m)', b);
%      else
%         error('Unrecognized mode %s', mode);
%      end
%   end
