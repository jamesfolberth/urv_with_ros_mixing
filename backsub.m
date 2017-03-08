function [x] = backsub(U,b)
%JMF: Do not use this for timing experiments!
   % column oriented, from G+VL 3.1.4
   [n,~] = size(b);
   x = b;
   for i=n:-1:2
      x(i,:) = x(i,:)/U(i,i);
      x(1:i-1,:) = x(1:i-1,:) - U(1:i-1,i)*x(i,:);
   end
   x(1,:) = x(1,:)/U(1,1);
end

