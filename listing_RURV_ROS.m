function [U,R,V,Vt] = RURV_ROS(A, n_its)
% RURV_ROS  RURV with ROS mixing for real matrices

   [m,n] = size(A);
   D_diags = sign(rand(n,n_its)-0.5); % diagonals of D_i; i.i.d. uniform +- 1
   
   % ROS mixing Ahat = A*V' using "transpose trick"
   Ahat = apply_V(A',D_diags)';

   % pre-sort
   nrms = sqrt(sum(Ahat.^2,1));
   [nrms,p] = sort(nrms, 2, 'descend');
   p_inv(p) = 1:numel(p);
   Ahat = Ahat(:,p);

   % unpivoted QR factorization
   [U,R] = qr(Ahat,0);
   
   % Return function handles to apply V and V^T on the left
   V = @(A) apply_V(A,D_diags,p);
   Vt = @(A) apply_Vt(A,D_diags,p_inv);
end

function [Ahat] = apply_V(A, D_diags, p)
% apply_V  Apply ROS mixing: V*A
   Ahat = A;
   for i=1:size(D_diags,2)
      Ahat = bsxfun(@times, Ahat, D_diags(:,i)); % Ahat = D_i*Ahat
      Ahat = dct(Ahat);                       % Ahat = F*Ahat;
   end

   if nargin == 3 % apply sorting
      Ahat = Ahat(p,:);
   end
end

function [Ahat] = apply_Vt(A, D_diags, p_inv)
% apply_Vt  Apply transpose ROS mixing: V^T*A
   if nargin == 3 % apply sorting
      Ahat = A(p_inv,:);
   else
      Ahat = A;
   end

   for i=size(D_diags,2):-1:1
      Ahat = idct(Ahat);                      % Ahat = F^T*Ahat;
      Ahat = bsxfun(@times, Ahat, D_diags(:,i)); % Ahat = D_i*Ahat
   end
end
