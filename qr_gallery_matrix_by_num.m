function [A] = qr_gallery_matrix_by_num(ind, n)
% Matrices from Table 2 of Demmel, Grigori, Gu, and Xiang's 
% "COMMUNICATION AVOIDING RANK REVEALING QR FACTORIZATION WITH COLUMN PIVOTING"

switch ind
case 1,
   A = qr_gallery('baart', n);
case 2,
   A = qr_gallery('break1', n);
case 3,
   A = qr_gallery('break9', n);
case 4,
   A = qr_gallery('deriv2', n);
case 5,
   A = qr_gallery('exponential', n);
case 6,
   A = qr_gallery('foxgood', n); % n must be divisible by 4
case 7,
   A = qr_gallery('gks', n);
case 8,
   A = qr_gallery('gravity', n);
case 9,
   A = qr_gallery('hc', n);
case 10,
   A = qr_gallery('heat', n);
case 11,
   A = qr_gallery('phillips', n);
case 12
   A = 2*rand(n)-1;
case 13
   A = 2*rand(n)-1; A = bsxfun(@times, A, 10*eps(1).^(((1:n)./n).'));
case 14,
   A = qr_gallery('shaw', n);
case 15,
   A = qr_gallery('spikes', n);
case 16,
   A = qr_gallery('stewart2', n);
case 17,
   A = qr_gallery('ursell', n);
case 18,
   A = qr_gallery('wing', n);
case 19,
   A = qr_gallery('kahan', n);
case 20,
   A = qr_gallery('devil', n);
otherwise
   error('Unknown index %d', ind);
end

end
