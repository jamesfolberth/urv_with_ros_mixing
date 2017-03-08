function [s] = dgejsv(A)
% DGEJSV  Call's LAPACK's dgejsv and returns the singular values
% This function uses Timothy Toolan's lapack mex file from 
% https://www.mathworks.com/matlabcentral/fileexchange/16777-lapack

   proto = 'DGEJSV(h,h,h,h,h,h,i,i,d,i,D,d,i,d,i,D,i,i,I)';
   [m,n] = size(A);
   
   if m == 1 && n == 1
      s = svd(A); % lapack segfaults if m==n==1
      return
   end

   A_copy = A;
   s = zeros(min(m,n),1);
   wrk = zeros(25*max(m,n),1);
   iwrk = zeros(m+3*n);
   [s,wrk,info] = lapack(proto, 'F', 'W', 'N', 'N', 'T', 'N', m, n, A_copy, m, s, zeros(m*n,1), m, [], 0, wrk, numel(wrk), iwrk, 0);
   s = wrk(1)/wrk(2)*s;

   if info > 0
      warning(strcat('DGEJSV did not converge in the given number of sweeps.',...
              '  The output singular values may be inaccurate.'));
   end
end


