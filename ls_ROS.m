function [x] = ls_ROS(A,b,n_its,underdet_mode,mix_fftw_opt)
% LS_ROS  This is our driver to solve LS problems with RURV_ROS/RVLU_ROS
   if nargin < 3, n_its = 1; end
   if nargin < 4, underdet_mode = 'basic'; end
   if nargin < 5, mix_fftw_opt = struct(); end
   
   [m,n] = size(A);
   do_VLU = 0;
   if m < n && strcmp(underdet_mode, 'minnorm')
      do_VLU = 1;
   end
   t = tic();

   opt_mix = mix_fftw_opt;
   opt_unmix = opt_mix;
   
   if do_VLU % want to mix A_mix^T = A^T*V^T <==> A_mix = V*A -> A=V'*L*U
      % to find the minnorm solution, we do a QR of A'
      % so A_mix = V*A, so the columns of A' are mixed
      % then sort A_mix' <- A_mix'*Pi' <==> A_mix <- Pi*A_mix
      % A = V'*Pi'*L*U 
      % A*x = b -> V'*Pi'*L*U*x = b -> L*U*x = Pi'*V'*b
      opt_mix.transpose = 0;
      opt_mix.side = 'l';

      [A_mix, darr] = mix_fftw(A,n_its,opt_mix);
      mix_time = toc(t);
      
      t = tic();
      [~,p] = sort(sum(A_mix.*conj(A_mix),2), 'descend'); % post-mixing naive pre-order
      A_mix = A_mix(p,:); 
      
      b_mix = mix_fftw(b, darr, opt_mix);
      b_mix = b_mix(p,:);
      unmix_time = toc(t);
      
      t = tic();
      x = qr_ls(A_mix, b_mix, underdet_mode);
      solve_time = toc(t);

   else % do a URV
      opt_mix.transpose = 1;
      opt_mix.side = 'r';
      opt_unmix.transpose = 1;
      opt_unmix.side = 'l';
   
      [A_mix, darr] = mix_fftw(A,n_its,opt_mix);
      mix_time = toc(t);
      
      t = tic();
      [~,p] = sort(sum(A_mix.*conj(A_mix),1), 'descend'); % post-mixing naive pre-order
	   pinv(p) = 1:numel(p);
      A_mix = A_mix(:,p); 
      
      y = qr_ls(A_mix, b, underdet_mode);
      solve_time = toc(t);
      
      t = tic();
      x = mix_fftw(y(pinv,:), darr, opt_unmix);
      unmix_time = toc(t);

   end
   
   verbose=1;
   if verbose
      fprintf(strcat('ls_ROS:\n  mix time = %3.2e\n  solve time = %3.2e\n',...
               '  unmix time = %3.2e\n  total time = %3.2e\n'),...
              mix_time, solve_time, unmix_time, mix_time+solve_time+unmix_time);
   end
  
end

