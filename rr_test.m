function [] = rr_test()

if ~add_rrqr_comparison(), return; end % for DGEQPX
if ~add_lapack(), return; end % for interface to LAPACK (DGEJSV)

%test_rr_conditions()
test_rr_condition_scaling()
%test_rank_approx()

end

function [tests] = ratios1(s,k,R11,R12,R22)
   % We use DGEJSV, which uses preconditioned Jacobi SVD, which should do better
   % on A = D1*C*D2, where D1,D2 diagonal and ill-cond, C well-cond
   sR11 = dgejsv(R11);
   sk_rat = max(1,s(k)/sR11(end));
   sR22 = dgejsv(R22);
   skp1_rat = max(1,sR22(1)/s(k+1));
   cond_R11 = cond(R11);
   norm_R11R12 = norm(R11\R12);
   tests = [sk_rat; skp1_rat; cond_R11; norm_R11R12];
end

function [tests] = ratios2(s,k,R11,R12,R22)
   % We use DGEJSV, which uses preconditioned Jacobi SVD, which should do better
   % on A = D1*C*D2, where D1,D2 diagonal and ill-cond, C well-cond
   sR11 = dgejsv(R11);
   sk_rat = max(1,max(s(1:k)./sR11(1:end)));
   sR22 = dgejsv(R22);
   skp1_rat = max(1,max(sR22(1:end)./s(k+1:end)));
   cond_R11 = cond(R11);
   norm_R11R12 = norm(R11\R12);
   tests = [sk_rat; skp1_rat; cond_R11; norm_R11R12];
end

function [A,s] = gen_mat(m,n,k)
   %[A,s] = qr_gallery('rankk', [m n], k);

   A = qr_gallery('kahan', [m n], acos(0.1));
   %A = qr_gallery('kahan_rankk', [m n], k, acos(0.1));
   %s = svd(A); % this can be quite inaccurate
   % We use DGEJSV, which uses preconditioned Jacobi SVD, which should do better
   % on A = D1*C*D2, where D1,D2 diagonal and ill-cond, C well-cond
   s = dgejsv(A);
end



%TODO it might be better to have the y-axis start at 1 and go to the max value across all 3 subplots
function [] = test_rr_conditions()
   %rng(271828);
   
   ratios = @(s,k,R11,R12,R22) ratios1(s,k,R11,R12,R22);
   
   m = 500; n = 500;
   %k = n-1;
   k = round(0.5*n);
   %[A,s] = qr_gallery('rankk', [m n], k);
   [A,s] = gen_mat(m,n,k);
   
   n_samples = 25;
   new_mat_per_sample = 0;
   names = {};
   rr_tests = zeros(4, numel(names), n_samples);

   for i=1:n_samples
      if new_mat_per_sample
         [A,s] = gen_mat(m,n,k);
      end
      cnt = 0;
      
      % QR
      cnt = cnt+1;
      names{cnt} = 'QR';
      [Q,R] = qr(A,0);
      R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
      rr_tests(:,cnt,i) = ratios(s,k,R11,R12,R22); 
      
      % QRCP
      % This is NOT a RRQR
      cnt = cnt+1;
      names{cnt} = 'QRCP';
      [Q,R,e] = qr(A,0);
      R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
      rr_tests(:,cnt,i) = ratios(s,k,R11,R12,R22); 
      
      % QR - From PGM
      % This is NOT a RRQR
      cnt = cnt+1;
      names{cnt} = 'QR - PGM';
      [Q,R,e] = qr_pgm_block(A);
      R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
      rr_tests(:,cnt,i) = ratios(s,k,R11,R12,R22); 

      % DGEQPX - From ACM algorithm 782
      % This is a RRQR
      cnt = cnt+1;
      names{cnt} = 'DGEQPX';
      [~,R,~] = qrxp(A);
      %[Q,R,e] = qrxp(A,0,0,0,size(A,1),2,s);
      R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
      rr_tests(:,cnt,i) = ratios(s,k,R11,R12,R22); 

      
      % URV Haar
      % This is WHP a RR URV
      cnt = cnt+1;
      names{cnt} = 'URV - Haar';
      [V,~] = qr(randn(n),0);
      [U,R] = qr(A*V',0);
      R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
      rr_tests(:,cnt,i) = ratios(s,k,R11,R12,R22); 

      % URV ROS
      cnt = cnt+1;
      names{cnt} = 'URV - ROS';
      mix_opt = mix_fftw();
      mix_opt.transpose=1;
      mix_opt.side='r';
      [Ahat,darr] = mix_fftw(A,10,mix_opt);
      [~,p] = sort(sum(Ahat.*conj(Ahat),1), 'descend'); % post-mixing naive pre-order
      Ahat = Ahat(:,p); 
      [U,R] = qr(Ahat,0);
      R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
      rr_tests(:,cnt,i) = ratios(s,k,R11,R12,R22); 
         
      if 0
         print_stuff = @(name, tests) fprintf(1,...
            ['%s\n',...
             '  s(k)/s_min(R11)    = %1.3e\n',...
             '  s_max(R22)/s(k+1)  = %1.3e\n',...
             '  norm(inv(R11)*R12) = %1.3e\n'], name, tests(1), tests(2), tests(4));
         for n=1:numel(names)
            print_stuff(names{n}, rr_tests(:,n,i));
         end
      end
   end

   rr_tests_mean = mean(rr_tests,3);
   rr_tests_std = std(rr_tests,[],3);

   print_stuff = @(name, tm, ts) fprintf(1,...
      ['%s\n',...
       '  s(k)/s_min(R11)    = %1.3e +- %1.3e\n',...
       '  s_max(R22)/s(k+1)  = %1.3e +- %1.3e\n',...
       '  norm(inv(R11)*R12) = %1.3e +- %1.3e\n'], name, tm(1), ts(1), tm(2), ts(2), tm(4), ts(4));
   for ni=1:numel(names)
      print_stuff(names{ni}, rr_tests_mean(:,ni), rr_tests_std(:,ni));
   end

   
   figure(1); clf;
   % https://www.mathworks.com/matlabcentral/newsreader/view_thread/171543 
   % https://www.mathworks.com/matlabcentral/fileexchange/45049-bland-altman-and-correlation-plot/content/BlandAltman/suptitle.m
   suptitle(sprintf('Samples of Rank-Revealing Conditions - m=%d, n=%d, k=%d',m,n,k));
   subplot(3,1,1);
   boxplot(permute(squeeze(rr_tests(1,:,:)),[2 1]), names);
   set(gca, 'yscale', 'log');
   m = ceil(log10(max(max(rr_tests(1,:,:)))));
   ytick = 10.^linspace(1,10*m,10*m);
   set(gca, 'YTick', ytick);
   %ylabel('\sigma(k)/\sigma_{min}(R_{11})');
   ylabel('\sigma(k)/\sigma_{min}(R_{11})');

   subplot(3,1,2);
   boxplot(permute(squeeze(rr_tests(2,:,:)),[2 1]), names);
   set(gca, 'yscale', 'log');
   m = ceil(log10(max(max(rr_tests(2,:,:)))));
   ytick = 10.^linspace(1,10*m,10*m);
   set(gca, 'YTick', ytick);
   %ylabel('s_max(R22)/s(k+1)', 'interpreter', 'none');
   ylabel('\sigma_{max}(R_{22})/\sigma(k+1)');

   subplot(3,1,3);
   boxplot(permute(squeeze(rr_tests(3,:,:)),[2 1]), names);
   set(gca, 'yscale', 'log');
   m = ceil(log10(max(max(rr_tests(3,:,:)))));
   ytick = 10.^linspace(1,10*m,10*m);
   set(gca, 'YTick', ytick);
   %ylabel('norm(inv(R11)*R12)', 'interpreter', 'none');
   ylabel('norm(R_{11}^{-1}R_{12})');

end

% how do max/avg/min values scale for kahan_rankk as n increases?
function [] = test_rr_condition_scaling()
   %rng(271828);
   
   ratios = @(s,k,R11,R12,R22) ratios2(s,k,R11,R12,R22);
   
   %mvec = round(logspace(1,log10(400),50));
   mvec = round(logspace(1,3,50));
   %mvec = round(logspace(1,4,50));
   %mvec = round(linspace(10,250,5));
   %kvec = round(0.5*mvec);
   kvec = mvec-1;
   
   n_samples = 5;
   new_mat_per_sample = 1;
   %names = {'QR', 'QRCP', 'QR - PGM', 'DGEQPX', 'URV - Haar', 'URV - ROS'};
   names = {'QRCP', 'HQRRP', 'DGEQPX', 'RURV\_Haar', 'RURV\_ROS'};
   rr_tests = zeros(4, numel(names), numel(mvec), n_samples);
      
   for i=1:numel(mvec)
      fprintf('\ri = %d',i);
      m = mvec(i); n = m;
      k = kvec(i);
      [A,s] = gen_mat(m,n,k);

      for si=1:n_samples
         if new_mat_per_sample && si > 1
            [A,s] = gen_mat(m,n,k);
         end
         
         % QR
         I = find(strcmp('QR', names));
         if ~isempty(I)
            ind = I(1);
            [Q,R] = qr(A,0);
            R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
            rr_tests(:,ind,i,si) = ratios(s,k,R11,R12,R22); 
         end
         
         % QRCP
         % This is NOT a RRQR
         I = find(strcmp('QRCP', names));
         if ~isempty(I)
            if n > 2000
               rr_tests(:,ind,i,si) = [1;1;1;1];
            else
               ind = I(1);
               [Q,R,e] = qr(A,0);
               R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
               rr_tests(:,ind,i,si) = ratios(s,k,R11,R12,R22);
            end
         end
         
         % QR - From PGM
         % This is NOT a RRQR
         I = find(strcmp('HQRRP', names));
         if ~isempty(I)
            ind = I(1);
            [Q,R,e] = qr_pgm_block(A);
            R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
            rr_tests(:,ind,i,si) = ratios(s,k,R11,R12,R22);
         end

         % DGEQPX - From ACM algorithm 782
         % This is a RRQR
         I = find(strcmp('DGEQPX', names));
         if ~isempty(I)
            ind = I(1);
            [~,R,~] = qrxp(A);
            %[Q,R,e] = qrxp(A,0,0,0,size(A,1),2,s);
            R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
            rr_tests(:,ind,i,si) = ratios(s,k,R11,R12,R22);
         end

         
         % URV Haar
         % This is WHP a RR URV
         I = find(strcmp('RURV\_Haar', names));
         if ~isempty(I)
            ind = I(1);
            [V,~] = qr(randn(n),0);
            [U,R] = qr(A*V',0);
            R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
            rr_tests(:,ind,i,si) = ratios(s,k,R11,R12,R22);
         end

         % URV ROS
         I = find(strcmp('RURV\_ROS', names));
         if ~isempty(I)
            ind = I(1);
            mix_opt = mix_fftw();
            mix_opt.transpose=1;
            mix_opt.side='r';
            [Ahat,darr] = mix_fftw(A,10,mix_opt);
            [~,p] = sort(sum(Ahat.*conj(Ahat),1), 'descend'); % post-mixing naive pre-order
            Ahat = Ahat(:,p); 
            [U,R] = qr(Ahat,0);
            R11 = R(1:k,1:k); R12 = R(1:k,k+1:n); R22 = R(k+1:m,k+1:n);
            rr_tests(:,ind,i,si) = ratios(s,k,R11,R12,R22);
         end
      end
   end
   fprintf('\r\n');

   rr_tests_mean = mean(rr_tests, 4);
   rr_tests_std = std(rr_tests, [], 4);
   rr_tests_max = max(rr_tests, [], 4);


   figure(1); clf;
   % https://www.mathworks.com/matlabcentral/newsreader/view_thread/171543 
   % https://www.mathworks.com/matlabcentral/fileexchange/45049-bland-altman-and-correlation-plot/content/BlandAltman/suptitle.m
   suptitle('Average of Sampled Rank-Revealing Conditions', 'fontsize', 30); %TODO not bold; not the same as title(..., 'fontsize', 20)
   
   subplot(3,1,1);
   holdMarker('reset');
   hold on;
   for i=1:numel(names)
      marker = holdMarker();
      plot(mvec, squeeze(rr_tests_max(1,i,:)), strcat('-',marker));
   end
   hold off;
   set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', 20);
   m = ceil(log10(max(max(max(rr_tests(1,:,:,:))))));
   %ytick = 10.^linspace(1,10*m,10*m);
   %set(gca, 'YTick', ytick);
   ylabel('max(\sigma_{i}(A)/\sigma_{i}(R_{11}))');
   legend(gca, deal(names), 'location', 'NorthWest');
   %axis([10 400 1 1e20])

   subplot(3,1,2);
   holdMarker('reset');
   hold on;
   for i=1:numel(names)
      marker = holdMarker();
      plot(mvec, squeeze(rr_tests_max(2,i,:)), strcat('-',marker));
   end
   hold off;
   set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', 20);
   m = ceil(log10(max(max(max(rr_tests(2,:,:,:))))));
   %ytick = 10.^linspace(1,10*m,10*m);
   %set(gca, 'YTick', ytick);
   ylabel('max(\sigma_{k+j}(A)/\sigma_{j}(R_{22}))');
   %axis([10 400 1 1e20])

   subplot(3,1,3);
   holdMarker('reset');
   hold on;
   for i=1:numel(names)
      marker = holdMarker();
      plot(mvec, squeeze(rr_tests_max(3,i,:)), strcat('-',marker));
   end
   hold off;
   set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', 20);
   m = ceil(log10(max(max(max(rr_tests(3,:,:,:))))));
   %ytick = 10.^linspace(1,10*m,10*m);
   %set(gca, 'YTick', ytick);
   xlabel('m');
   ylabel('norm(R_{11}^{-1}R_{12})');
   %axis([10 400 1 1e20])


   %subplot(3,1,2);
   %boxplot(permute(squeeze(rr_tests(2,:,:)),[2 1]), names);
   %set(gca, 'yscale', 'log');
   %m = ceil(log10(max(max(rr_tests(2,:,:)))));
   %ytick = 10.^linspace(1,10*m,10*m);
   %set(gca, 'YTick', ytick);
   %%ylabel('s_max(R22)/s(k+1)', 'interpreter', 'none');
   %ylabel('\sigma_{max}(R_{22})/\sigma(k+1)');

   %subplot(3,1,3);
   %boxplot(permute(squeeze(rr_tests(3,:,:)),[2 1]), names);
   %set(gca, 'yscale', 'log');
   %m = ceil(log10(max(max(rr_tests(3,:,:)))));
   %ytick = 10.^linspace(1,10*m,10*m);
   %set(gca, 'YTick', ytick);
   %%ylabel('norm(inv(R11)*R12)', 'interpreter', 'none');
   %ylabel('norm(R_{11}^{-1}R_{12})');



end

function [] = test_rank_approx()
   
   m = 1000;
   A = qr_gallery('kahan', m, acos(0.1));
   %A = qr_gallery('rankk', 100, 50);
   
   %S = load('sjsu_smdb/lp_bore3d.mat');
   %S = load('sjsu_smdb/NotreDame_yeast.mat');
   %S = load('sjsu_smdb/lock2232.mat'); % note that `svd` doesn't match S.Problem.svals!
   %A = full(S.Problem.A);

   %s = svd(A);
   s = dgejsv(A);
   %s = S.Problem.svals;
   
   names = {};
   kvec = floor(linspace(1, 80, 25)).';
   %kvec = floor(linspace(1, min(size(A))-1, 25)).';
   %kvec = floor(linspace(1, 500, 25)).';
   errors = zeros(numel(kvec),2,1);
   
   xax = [];
   %xax = [0 2232];
   yax = [];
   yax = [1e0 1e3];
   
   do_specnorm = 1;
   do_frobnorm = 1;
   
   %TODO should we do QR, RSVD here too?

   ind = 0;

   % QR
   ind = ind + 1;
   names{ind} = 'QR';
   [Q,R] = qr(A,0);
   for ki=1:numel(kvec)
      k = kvec(ki)
      err = A - Q(:,1:k)*R(1:k,:);
      if do_specnorm, errors(ki,1,ind) = norm(err); end
      if do_frobnorm, errors(ki,2,ind) = norm(err,'fro'); end
   end

   % QRCP
   ind = ind + 1;
   names{ind} = 'QRCP';
   [Q,R,e] = qr(A,0);
   for ki=1:numel(kvec)
      k = kvec(ki)
      err = A(:,e) - Q(:,1:k)*R(1:k,:);
      if do_specnorm, errors(ki,1,ind) = norm(err); end
      if do_frobnorm, errors(ki,2,ind) = norm(err,'fro'); end
   end

   % QR - PGM
   ind = ind + 1;
   names{ind} = 'QRPGM';
   [Q,R,e] = qr_pgm_block(A);
   for ki=1:numel(kvec)
      k = kvec(ki)
      err = A(:,e) - Q(:,1:k)*R(1:k,:);
      if do_specnorm, errors(ki,1,ind) = norm(err); end
      if do_frobnorm, errors(ki,2,ind) = norm(err,'fro'); end
   end

   % DGEQPX
   ind = ind + 1;
   names{ind} = 'DGEQPX';
   %[Q,R,e] = qr_pgm_block(A);
   [Q,R,P] = qrxp(A,0,0,0,size(A,1),2);
   for ki=1:numel(kvec)
      k = kvec(ki)
      err = A*P - Q(:,1:k)*R(1:k,:);
      if do_specnorm, errors(ki,1,ind) = norm(err); end
      if do_frobnorm, errors(ki,2,ind) = norm(err,'fro'); end
   end
   
   % QR - Haar
   ind = ind + 1;
   names{ind} = 'URVHaar';
   [V,~] = qr(randn(size(A,2)),0);
   [U,R] = qr(A*V',0);
   for ki=1:numel(kvec)
      k = kvec(ki)
      err = A - U(:,1:k)*R(1:k,:)*V;
      if do_specnorm, errors(ki,1,ind) = norm(err); end
      if do_frobnorm, errors(ki,2,ind) = norm(err,'fro'); end
   end

   % QR - ROS
   ind = ind + 1;
   names{ind} = 'URVROS';
   mix_opt = mix_fftw();
   mix_opt.transpose=1;
   mix_opt.side = 'r';
   unmix_opt = mix_opt;
   unmix_opt.side = 'l';

   [Ahat,darr] = mix_fftw(A,1,mix_opt);
   [U,R] = qr(Ahat,0);
   for ki=1:numel(kvec)
      k = kvec(ki)
      err = A - mix_fftw(U(:,1:k)*R(1:k,:), darr, unmix_opt);
      if do_specnorm, errors(ki,1,ind) = norm(err); end
      if do_frobnorm, errors(ki,2,ind) = norm(err,'fro'); end
   end

   figure(1); clf;
   if isempty(xax)
      xmin = min(kvec); xmax = max(kvec);
   else
      xmin = xax(1); xmax = xax(2);
   end

   %if do_specnorm && do_frobnorm, subplot(1,2,1); end
   if do_specnorm
      hold on
      plot(kvec, s(kvec+1), '-k');
      for i=1:numel(names)
         plot(kvec, errors(:,1,i), '-o');
      end
      hold off
      set(gca, 'yscale', 'log');
      %if ~do_frobnorm
         legend(deal(horzcat({'DGEJSV'}, names)), 'Location', 'NorthEast');
      %end
      ylabel('norm(A-A_k)');
      xlabel('Reconstruction Rank - k');
      if ~isempty(yax)
         axis([xmin xmax yax(1) yax(2)]);
      end
   end


   %if do_specnorm && do_frobnorm, subplot(1,2,2); end
   figure(2); clf;
   if do_frobnorm
      hold on
      frob_errs = sqrt(cumsum(s(end:-1:1).^2)); frob_errs = frob_errs(end:-1:1);
      plot(kvec, frob_errs(kvec+1), '-k');
      for i=1:numel(names)
         plot(kvec, errors(:,2,i), '-o');
      end
      hold off
      set(gca, 'yscale', 'log');
      legend(deal(horzcat({'DGEJSV'}, names)), 'Location', 'NorthEast');
      ylabel('norm(A-A_k,''fro'')');
      xlabel('Reconstruction Rank - k');
      if ~isempty(yax)
         axis([xmin xmax yax(1) yax(2)]);
      end
   end

end


function [Q,R,p] = qr_pgm_block(A)
	% {{{ blocked QR (normal pivot) from Nathan/Gunnar
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Given an n x n matrix A, this function produces a factorization
	%    A(:,J) = Q * R
	% where 
	%    J is a permutation vector of length n,
	%    Q is an n x n orthonormal matrix,
	%    R is an n x n upper triangular matrix.
	%
	% It operates in blocked fashion, processing "b" columns at a time.
	%
	% The input parameter "p" indicates how much over-sampling we do.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	function [Q,R,J] = LOCAL_qr(A,b,p)
	
	n     = size(A,2);
	nstep = ceil(n/b);
	Q     = eye(n);
	J     = 1:n;
	
	%%% We will create an upper triangular matrix "R" by applying a sequence
	%%% of ON maps to the matrix "A". For now, we simply copy A to R, and we 
	%%% will then perform the various transforms (overwriting R at each step).
	R = A;
	
	%%% Process all blocks, except the last.
	%%% At each step, we partition the index vector as follows:
	%%%    (1:n) = [I1, I2, I3]
	%%% where "I2" is the block currently being processed.
	for j = 1:nstep
	
	  %%% Construct the index vectors that partition the matrix.  
	  I1 = 1:((j-1)*b);
	  I2 = (j-1)*b + (1:min(b,n-b*(j-1)));
	  I3 = (j*b+1):n;
	  
	  %%% Find b good pivot columns from [I2 I3] using the randomized sampling
	  %%% procedure and move them to the I2 column.
	  %%% (We don't do this at the last step, when I3 is empty.)
	  if (j < nstep)
	    %%% Compute the "sampling matrix" Y.
	    G = randn(b+p,n-(j-1)*b);
	    Y = G * R([I2,I3],[I2,I3]);
	    %%% Execute full pivoted QR on the small matrix Y (only b steps),
	    [~,~,Jloc] = pivoted_QR(Y,b);
	    %%% Permute the columns in the [I2,I3] block as dictated by Jloc:
	    I23          = [I2 I3];
	    R(:,[I2,I3]) = R(:,I23(Jloc));
	    J([I2,I3])   = J(I23(Jloc));  
	  end
	  
	  %%% Next execute a full QR on the middle block.
	  %%% Horrible inefficiency here since Qloc is a product of b Householder
	  %%% reflectors!!!
	  [Qloc,Rloc,Jloc] = qr(R([I2,I3],I2),'vector');
	  R([I2,I3],I2)    = Rloc;
	  R(I1,I2)         = R(I1,I2(Jloc));
	  R([I2,I3],I3)    = Qloc'*R([I2,I3],I3);  % This is a rank-b update!!!
	  Q(:,[I2,I3])     = Q(:,[I2,I3])*Qloc;    % This is a rank-b update!!!
	  J(I2)            = J(I2(Jloc));
	  
	end
	
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function performs classical column pivoted QR.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	function [Q,R,ind] = pivoted_QR(A,k)
	
	% This function orthogonalizes the COLUMNS of A
	% It uses a modified Gram-Schmidt with column pivoting
	
	m = size(A,1);
	n = size(A,2);
	
	R = zeros(min(m,n),n);
	Q = A;
	ind = 1:n;
	
	for j = 1:k
	    [~, j_max] = max(sum(Q(:,j:n).*Q(:,j:n)));
	    j_max          = j_max + j - 1;
	    Q(:,[j,j_max]) = Q(:,[j_max,j]);
	    R(:,[j,j_max]) = R(:,[j_max,j]);
	    ind([j,j_max]) = ind([j_max,j]);
	    r_jj   = norm(Q(:,j));
	    Q(:,j) = Q(:,j)/r_jj;
	    Q(:,j) = Q(:,j) - Q(:,1:(j-1))*(Q(:,1:(j-1))'*Q(:,j));
	    Q(:,j) = Q(:,j)/norm(Q(:,j));
	    R(j,j) = r_jj;
	    rr     = Q(:,j)'*Q(:,(j+1):n);
	    R(j,(j+1):n) = rr;
	    Q(:,(j+1):n) = Q(:,(j+1):n) - Q(:,j)*rr;
	end
	
	Q = Q(:, 1:min(m,n));
	
	end

	% }}}
	[Q,R,p] = LOCAL_qr(A,20,10);
end


