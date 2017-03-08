function [] = ls_test(times_file, times_file2)
% LS_TEST  Various experiments with least-squares
%
% We set up parameters of the experiments by modifying appropriate driver functions.
%
% For ls_timing(), we can call ls_test with one or two arguments to plot timing results.
% The syntax and behavior is similar to qr_timing.m

if nargin == 1
   ls_timing_plot(times_file);
   return
elseif nargin == 2
   ls_timing_plot2(times_file, times_file2);
   return
end

% Check for BLENDENPIK in the path
if ~exist('blendenpik')
   error(sprintf(['Please install, configure, and add BLENDENPIK to your path.\n'...
          'You can download BLENDENPIK from\n'...
          'https://www.mathworks.com/matlabcentral/fileexchange/25241-blendenpik']));
end

%simple_test();
%ls_timing_single()
ls_timing()
%ls_correlated_basic();

end


function [] = simple_test()
   rng(271828);
   n = 500;
   m = 10*n;
   m=50; n=10;
   %m=5e4; n=5e2;
   %m=1000; n=10*m;
   %m=1e2; n=5e4; % 'basic' beats BLENDENPIK by a bit
   k = 25;
   
   % From Stephen's initial experiments
   A = randn(m,n) + exp(10*rand(m,n));
   % If we do this, then more obvious difference///
   A = bsxfun( @times, A, exp(2*rand(1,n)) );
   %A = qr_gallery('kahan', n);
   %A = qr_gallery('extkahan', n);
   %[A,s] = qr_gallery('devil', n, 50);
   %A = qr_gallery('stewart'); % this gives some large deviations
   %[A,s] = qr_gallery('rankk', [m n], 25);
   A = qr_gallery('condk', [m n], 1e6);
   
   b = randn(m,1);
   
   t = tic();
   x_ls = A\b; % if underdetermined, \ will not give the minimum norm solution.
   %x_ls = pinv(A)*b; % will produce minimum norm solution
   bs_time = toc(t);

   %x_ls = lsqr(A,b,1e-14,100);
   Aerr = @(x) norm(A*(x-x_ls));
   resid = @(x) norm(A*x-b);

   mix_fftw_opt = mix_fftw();
	mix_fftw_opt.plan_rigor = 1;
	mix_fftw_opt.wisdom_file = 'wisdom/wisdom.fftw';
	mix_fftw_opt.nthreads = 2;
	mix_fftw_opt.block_size = 100;
   mix_fftw_opt.use_darr = 0;%XXX temp test

   t = tic();
   x_ls_ROS = ls_ROS(A,b,1,'basic',mix_fftw_opt);
   dct_time = toc(t);
   
   fprintf('\\ time       = %1.3e sec\n', bs_time);
   fprintf('Aerr (\\)     = %1.3e\n', Aerr(x_ls));
   fprintf('resid (\\)    = %1.3e\n', resid(x_ls));
   fprintf('norm (\\)     = %1.3e\n', norm(x_ls));

   fprintf('DCT time     = %1.3e sec\n', dct_time);
   fprintf('Aerr (DCT)   = %1.3e\n', Aerr(x_ls_ROS));
   fprintf('resid (DCT)  = %1.3e\n', resid(x_ls_ROS));
   fprintf('norm (DCT)   = %1.3e\n', norm(x_ls_ROS));
   
   % BLENDENPIK appears to be getting very close to minimum norm solution
   t = tic();
   x_ls_bpik = ls_blendenpik(A,b);
   bpik_time = toc(t);

   fprintf('BPIK time    = %1.3e sec\n', bpik_time);
   fprintf('Aerr (BPIK)  = %1.3e\n', Aerr(x_ls_bpik));
   fprintf('resid (BPIK) = %1.3e\n', resid(x_ls_bpik));
   fprintf('norm (BPIK)  = %1.3e\n', norm(x_ls_bpik));

   %x_ls_pgm = ls_pgm(A,b);
   %fprintf('Aerr (PGM)   = %1.3e\n', Aerr(x_ls_pgm));
   %fprintf('resid (PGM)  = %1.3e\n', resid(x_ls_pgm));
   %fprintf('norm (PGM)   = %1.3e\n', norm(x_ls_pgm));

end


function [] = ls_timing_single()
   
   m = 100;
   aspect = 0.1; % height/width (overdet is > 1; underdet is < 1);
   n_its = 1;
   n_samples = 10;
   
   n = round(m/aspect);
   A = qr_gallery('condk', [m n], 1e6);
   %b = randn(m,1)
   b = A*randn(n,1);

   mix_fftw_opt = mix_fftw();
	mix_fftw_opt.plan_rigor = 1;
	mix_fftw_opt.wisdom_file = '/scratch/wisdom/wisdom.fftw';
	mix_fftw_opt.nthreads = 2;
	mix_fftw_opt.block_size = 100;

   funcs = {'A\b', 'BLENDENPIK', 'ls_ROS'};
   times = zeros(n_samples, numel(funcs));

   for i=1:numel(funcs)
      func = funcs{i};
         
      if strcmp(func, 'A\b')
         x = A\b;
         for s=1:n_samples
            dt = tic();
            x = A\b;
            times(s,i) = toc(dt);
         end

      elseif strcmp(func, 'BLENDENPIK')
         x = blendenpik(A,b);
         for s=1:n_samples
            dt = tic();
            x = blendenpik(A,b);
            times(s,i) = toc(dt);
         end

      elseif strcmp(func, 'ls_ROS')
         x = ls_ROS(A,b,n_its,'basic',mix_fftw_opt);
         for s=1:n_samples
            dt = tic();
            x = ls_ROS(A,b,n_its,'basic',mix_fftw_opt);
            times(s,i) = toc(dt);
         end
      end
   end
   
   avg = squeeze(mean(times,1));
   sdev = squeeze(std(times,0,1));
   
   fprintf('\n');
   fprintf('Size of A is [%d %d]\n', m, n);
   for i=1:numel(funcs)
      fprintf('% 10s: %1.2e +- %1.2e (s)\n', funcs{i}, avg(i), sdev(i));
   end
      
end

function ls_timing()

   now = datetime('now');
   times_file = sprintf('ls_timing_times_%d%02d%02d_%02d%02d%02d.mat', now.Year, now.Month, now.Day,...
      now.Hour, now.Minute, round(now.Second));
   %maxNumCompThreads(1); %XXX
   n_threads = maxNumCompThreads();

   mvec = round(linspace(100,500,5));
  
   %mvec = round(linspace(100,5000,20)); % 1/20
   mvec = round(linspace(100,5000,20)); % 1/2
   %mvec = round(linspace(100,10000,20)); % 2
   %mvec = round(linspace(100,35000,20)); % 2
   %mvec = round(linspace(100,50000,20)); % 20
   %mvec = round(linspace(100,125000,25)); % aspect=20 and we're out of memory at 133333x6667
   aspect = 1/20; % m/n (overdet is > 1; underdet is < 1);
   n_its = 1;
   n_samples = 3;

   mix_fftw_opt = mix_fftw();
	mix_fftw_opt.plan_rigor = 4
	mix_fftw_opt.wisdom_file = '/scratch/wisdom/wisdom.fftw';
	mix_fftw_opt.nthreads = n_threads;
	mix_fftw_opt.block_size = 250;

   %funcs = {'A\b', 'BLENDENPIK', 'ls_ROS (basic)', 'ls_ROS (minnorm)'};
   if aspect < 1
      funcs = {'BLENDENPIK', 'ls_ROS (basic)', 'ls_ROS (minnorm)'};
      %funcs = {'ls_ROS (basic)', 'ls_ROS (minnorm)'}; % thread timing
      %funcs = {'A\b', 'BLENDENPIK', 'ls_ROS (basic)', 'ls_ROS (minnorm)'};
   else
      %funcs = {'BLENDENPIK', 'ls_ROS'}; 
      funcs = {'ls_ROS'}; % thread timing
   end
   times = zeros(n_samples, numel(funcs), numel(mvec));
   resids = zeros(n_samples, numel(funcs), numel(mvec));

   bpik_opts = struct('type', 'DCT');
  
   for mi=1:numel(mvec)
      fprintf('mi = %d of %d\n', mi, numel(mvec));
      m = mvec(mi);
      n = round(m/aspect);
      fprintf('[m n] = [%d %d]\n', m, n);
      A = qr_gallery('condk', [m n], 1e6);
      %b = randn(m,1)
      b = A*randn(n,1);

      for fi=1:numel(funcs)
         func = funcs{fi};
            
         if strcmp(func, 'A\b')
            x = A\b;
            for s=1:n_samples
               dt = tic();
               x = A\b;
               times(s,fi,mi) = toc(dt);
               resids(s,fi,mi) = norm(A*x-b);
               fprintf('% 20s: norm(A*x-b) = %1.3e\n', func, resids(s,fi,mi));
               fprintf('% 20s: norm(x)     = %1.3e\n\n', func, norm(x));
            end

         elseif strcmp(func, 'BLENDENPIK')
            x = blendenpik(A,b,bpik_opts);
            for s=1:n_samples
               dt = tic();
               x = blendenpik(A,b,bpik_opts);
               times(s,fi,mi) = toc(dt);
               resids(s,fi,mi) = norm(A*x-b);
               fprintf('% 20s: norm(A*x-b) = %1.3e\n', func, resids(s,fi,mi));
               fprintf('% 20s: norm(x)     = %1.3e\n\n', func, norm(x));
            end

         elseif strcmp(func, 'ls_ROS') || strcmp(func, 'ls_ROS (basic)')
            x = ls_ROS(A,b,n_its,'basic',mix_fftw_opt);
            for s=1:n_samples
               dt = tic();
               x = ls_ROS(A,b,n_its,'basic',mix_fftw_opt);
               times(s,fi,mi) = toc(dt);
               resids(s,fi,mi) = norm(A*x-b);
               fprintf('% 20s: norm(A*x-b) = %1.3e\n', func, resids(s,fi,mi));
               fprintf('% 20s: norm(x)     = %1.3e\n\n', func, norm(x));
            end

         elseif strcmp(func, 'ls_ROS (minnorm)')
            x = ls_ROS(A,b,n_its,'minnorm',mix_fftw_opt);
            for s=1:n_samples
               dt = tic();
               x = ls_ROS(A,b,n_its,'minnorm',mix_fftw_opt);
               times(s,fi,mi) = toc(dt);
               resids(s,fi,mi) = norm(A*x-b);
               fprintf('% 20s: norm(A*x-b) = %1.3e\n', func, resids(s,fi,mi));
               fprintf('% 20s: norm(x)     = %1.3e\n\n', func, norm(x));
            end

         end
      end
      %save(times_file, 'times', 'resids', 'n_samples', 'funcs', 'mvec', 'n_threads', 'aspect', 'n_its');
   end

   %times_file = sprintf('ls_timing_times_%d%02d%02d_%02d%02d%02d.mat', now.Year, now.Month, now.Day,...
   %   now.Hour, now.Minute, round(now.Second));
   %n_threads = maxNumCompThreads();
   save(times_file, 'times', 'resids', 'n_samples', 'funcs', 'mvec', 'n_threads', 'aspect', 'n_its');
   
   %ls_timing_plot(times_file);

end

function ls_timing_plot(times_file);
   
   S = load(times_file);
   times = S.times;
   resids = S.resids;
   n_samples = S.n_samples;
   funcs = S.funcs;
   funcs
   % override names
   funcs = {'BLENDENPIK', 'RURV_ROS (basic)', 'RVLU_ROS (minnorm)'}
   %funcs = {'BLENDENPIK', 'RURV_ROS'}
   mvec = S.mvec;
   n_threads = S.n_threads;
   aspect = S.aspect;
   
   plot_times = 1;
   plot_resids = 0;
   %TODO may want to plot avg-min, avg+max instead of avg+-std; this can be done on saved data
   
   if plot_times
      avg = squeeze(mean(times,1));
      sdev = squeeze(std(times,0,1));
      
      figure(1);
      clf;
      hold on
      for fi=1:numel(funcs)
         func = funcs{fi};
         %errorbar(mvec, avg(fi,:), sdev(fi,:),'-o');
         loglog(mvec, avg(fi,:), '-o');
      end
      hold off

      set(gca, 'yscale', 'log', 'xscale', 'log', 'fontsize', 40);
      if n_threads == 1
         title(sprintf('Average runtime for m/n = %1.2f for %d runs (%d thread)', aspect, n_samples, n_threads),'fontsize', 40);
      else
         title(sprintf('Average runtime for m/n = %1.2f for %d runs (%d threads)', aspect, n_samples, n_threads), 'fontsize', 40);
      end
      xlabel('m');
      ylabel('Average runtime (s)');
      legend(gca, deal(funcs), 'Location', 'NorthWest', 'interpreter', 'none');
      axis([1e2 5e3 1e-2 1e3]) % underdet
      %axis([1e2 5e4 1e-2 1e2]) % overdet
   end

   if plot_resids
      avg = squeeze(mean(resids,1));
      sdev = squeeze(std(resids,0,1));
      
      figure(2);
      clf;
      hold on
      for fi=1:numel(funcs)
         func = funcs{fi};
         errorbar(mvec, avg(fi,:), sdev(fi,:),'-o');
      end
      hold off

      set(gca, 'yscale', 'log');
      title(sprintf('Average residual (same A) for m/n = %1.2f for %d runs', aspect, n_samples));
      xlabel('m');
      ylabel('Average residual (using the same A)');
      legend(gca, deal(funcs), 'Location', 'NorthWest', 'interpreter', 'none');
   end
end

function ls_timing_plot2(times_file, times_file2);
   
   S = load(times_file);
   times = S.times;
   resids = S.resids;
   n_samples = S.n_samples;
   funcs = S.funcs;
   funcs
   mvec = S.mvec;
   n_threads = S.n_threads;
   aspect = S.aspect;

   S2 = load(times_file2);
   times2 = S2.times;
   resids2 = S2.resids;
   n_samples2 = S2.n_samples;
   funcs2 = S2.funcs;
   funcs
   mvec2 = S2.mvec;
   n_threads2 = S2.n_threads;
   aspect2 = S2.aspect;

   %assert(funcs == funcs2);
   %assert(mvec == mvec2);
   assert(aspect == aspect2);

   funcs = {'RURV_ROS (basic)', 'RVLU_ROS (minnorm)'}
   %funcs2 = funcs;
   %funcs = {'RURV_ROS'}
   funcs2 = funcs;
   
   %TODO may want to plot avg-min, avg+max instead of avg+-std; this can be done on saved data
   
   co = get(gca, 'colororder');
   avg = squeeze(mean(times,1));
   sdev = squeeze(std(times,0,1));

   avg2 = squeeze(mean(times2,1));
   sdev2 = squeeze(std(times2,0,1));
   
   plot_ratios = 1;
   

   if plot_ratios
      figure(1);
      clf;
      hold on
      if numel(funcs) > 1
         for fi=1:numel(funcs)
            func = funcs{fi};
            semilogx(mvec, avg(fi,:)./avg2(fi,:), '-o', 'color', co(fi,:));
         end
      else
         semilogx(mvec, avg2./avg, '-o', 'color', co(1,:));
      end
      hold off
         
      axis([1e2 5e3 0 10])
      %axis([1e2 5e4 0 10])
      
      labels = {};
      for i=1:numel(funcs)
         labels{i} = sprintf('%s - %d/%d threads', funcs{i}, n_threads2, n_threads)
      end

      legend(gca, deal(labels), 'Location', 'NorthWest', 'interpreter', 'none');
      title(sprintf('Speedup Factor for m/n = %1.2f for %d runs', aspect, n_samples), 'fontsize', 40);
      xlabel('m');
      ylabel('Speedup Factor');
      set(gca, 'xscale', 'log', 'fontsize', 40);
   
   else
      figure(1);
      clf;
      hold on
      if numel(funcs) > 1
         for fi=1:numel(funcs)
            func = funcs{fi};
            %errorbar(mvec, avg(fi,:), sdev(fi,:),'-o', 'color', co(:,fi));
            %errorbar(mvec, avg2(fi,:), sdev2(fi,:),'--x','color', co(:,fi));
            loglog(mvec, avg(fi,:), '-o', 'color', co(fi,:));
            loglog(mvec, avg2(fi,:),'--x', 'color', co(fi,:));
         end
      else
         loglog(mvec, avg, '-o', 'color', co(1,:));
         loglog(mvec, avg2,'--x', 'color', co(1,:));
      end
      hold off
         
      axis([1e2 5e3 5e-3 1e2])
      %axis([1e2 5e4 1e-2 1e3])

     for i=1:numel(funcs)
        if n_threads > 1
           funcs{i} = sprintf('%s - %d threads', funcs{i}, n_threads);
        else
           funcs{i} = sprintf('%s - %d thread', funcs{i}, n_threads);
        end
     end
     for i=1:numel(funcs2)
        if n_threads2 > 1
           funcs2{i} = sprintf('%s - %d threads', funcs2{i}, n_threads2);
        else
           funcs2{i} = sprintf('%s - %d thread', funcs2{i}, n_threads2);
        end
     end

      legend(gca, deal(horzcat(funcs, funcs2)), 'Location', 'NorthWest', 'interpreter', 'none');
      title(sprintf('Average runtime for m/n = %1.2f for %d runs', aspect, n_samples), 'fontsize', 40);
      xlabel('m');
      ylabel('Average runtime (s)');
      set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', 40);
   end

end


function [] = ls_correlated_basic()

   rng(271828);
   
   n_threads = maxNumCompThreads();
   
   n_samples = 10;

   n_its = 1;
   mix_fftw_opt = mix_fftw();
	mix_fftw_opt.plan_rigor = 4;
	mix_fftw_opt.wisdom_file = '/scratch/wisdom/wisdom.fftw';
	mix_fftw_opt.nthreads = n_threads;
	mix_fftw_opt.block_size = 250;

   m = 1000; n = 1500;
   %overlap_frac = .01;
   p = 10;
   e = 1e-4;
   %if round(overlap_frac) == overlap_frac
   %   A = repmat(randn(m,n), [1 ooverlap_frac]);
   %else
   %   ncol = max(0,round(overlap_frac*n))
   %   if overlap_frac >= 1, A = repmat(randn(m,n), [1 floor(overlap_frac)]); end
   %   A = randn(m,n);
   %   p = randperm(n);
   %   A = [A A(:,p(1:ncol))];
   %end
   %A = A + e*randn(size(A));
   %p = randperm(size(A,2));
   %%p = randperm(n);
   %num_repeated = sum(p(1:m) > n)
   %A = A(:,p);

   %A = randn(m,n);
   %perm = randperm(n,p);
   %A = [A A(:,perm)];
   %perm = randperm(n+p,n);
   %%num_repeated = sum(perm(1:m) > n)
   %A = A(:,perm) + e*randn(m,n);

   A = randn(m,n-p);
   %TODO adding this to randomly scale colums causes BLENDENPIK to have a
   %     resid on par with QR/QRCP.  RURVs are 2 orders better
   %A = bsxfun(@times, A, exp(10*rand(1,size(A,2))));
   perm = randperm(n-p,p);
   A = [A A(:,perm)]; % append 1st p columns
   perm = randperm(n);
   A = A(:,perm) + e*randn(m,n);


   %A = [1 0 0 0; 0 1 1 1; 0 0 1e-100 1];
   %[m,n] = size(A);

   cond(A)
   %svd(A)

   %[Q,R] = qr(A,0);
   %R
   %imagesc(abs(A))

   b = randn(m,1);
   X_true = ones(size(A,2),size(b,2)); % this is bogus

   resid = @(x) norm(A*x-b,'fro');
   errorf = @(x) norm(x-X_true,'fro');
   
   times_ = zeros(n_samples, 5);

   x_qr = qr_ls(A,b,'basic');
   for i=1:n_samples
      % QR without pivoting
      t = tic();
      x_qr = qr_ls(A,b,'basic');
      %[Q,R] = qr(A,0);
      %x_qr = backsub(R,Q'*b);
      %% basic solution
      %if size(R,1) < size(A,2), x_qr(size(R,1)+1:size(A,2),:) = 0; end
      time_qr = toc(t);
      norm_qr = norm(x_qr,'fro');
      resid_qr = resid(x_qr);
      error_qr = errorf(x_qr);
      times_(i,1) = time_qr;
   end   
      
   %x_qrcp = qrcp_ls(A,b,'basic');
   x_qrcp = qrcp_ls(A,b,'basicfull');
   for i=1:n_samples
      % \ (probably uses QRCP) or do our own QRCP
      t = tic();
      %x_qrcp = A\b;
      %QRCP
      %[Q,R,e] = qr(A,0);
      %x_qrcp = backsub(R, Q'*b);
      % basic solution
      %if size(R,1) < size(A,2), x_qrcp(size(R,1)+1:size(A,2),:) = 0; end
      %if size(A,1) < size(A,2), x_qrcp(size(A,1)+1:size(A,2),:) = 0; end
      %einv(e) = 1:numel(e);
      %x_qrcp = x_qrcp(einv,:);
      %if size(A,1) < size(A,2), x_qrcp(size(A,1)+1:size(A,2),:) = 0; end
      %x_qrcp = qrcp_ls(A,b,'basic');
      x_qrcp = qrcp_ls(A,b,'basicfull');
      %SVD
      %[U,S,V] = svd(A,'econ');
      %x_qrcp = V*((U'*b)./diag(S));
      time_qrcp = toc(t);
      norm_qrcp = norm(x_qrcp,'fro');
      resid_qrcp = resid(x_qrcp);
      error_qrcp = errorf(x_qrcp);
      times_(i,2) = time_qrcp;
   end
   
   x_bpik = ls_blendenpik(A,b);
   for i=1:n_samples
      % BLENDENPIK (for Kahan, just falls back on unpivoted QR in LAPACK xGELS)
      t = tic();
      % minnorm solution
      x_bpik = ls_blendenpik(A,b);
      % basic solution
      %x_bpik = ls_blendenpik(A(:,1:m),b);
      %if size(A,1) < size(A,2), x_bpik(size(A,1)+1:size(A,2),:) = 0; end
      time_bpik = toc(t);
      norm_bpik = norm(x_bpik,'fro');
      resid_bpik = resid(x_bpik);
      error_bpik = errorf(x_bpik);
      times_(i,3) = time_bpik;
   end
      
   for i=1:n_samples
      % Haar mixed URV
      t = tic();
      [V,~] = qr(randn(size(A,2)),0);
      Ah = A*V';
      if m < n
         [U,R] = qr(Ah(1:m,1:m),0);
      else
         [U,R] = qr(Ah,0);
      end
      x_haar = backsub(R,U'*b);
      % basic solution
      if size(R,1) < size(A,2), x_haar(size(R,1)+1:size(A,2),:) = 0; end
      x_haar = V'*x_haar;
      time_haar = toc(t);
      norm_haar = norm(x_haar,'fro');
      resid_haar = resid(x_haar);
      error_haar = errorf(x_haar);
      times_(i,4) = time_haar;
   end
   
   x_mix = ls_ROS(A,b,n_its,'basic',mix_fftw_opt); %TODO is mix_fftw not caching wisdom correctly?
   for i=1:n_samples
      % our mixed URV
      t = tic();
      x_mix = ls_ROS(A,b,n_its,'basic',mix_fftw_opt);
      time_mix = toc(t)
      norm_mix = norm(x_mix,'fro');
      resid_mix = resid(x_mix);
      error_mix = errorf(x_mix);
      times_(i,5) = time_mix;
   end 
   
   

   %print_stuff = @(name,time,r) fprintf('%s\n  time = %1.3e sec\n  resid = %1.3e\n',name,time,r);
   %print_stuff('QR', time_qr, resid_qr);
   %print_stuff('QRCP', time_qrcp, resid_qrcp);
   %print_stuff('BLENDENPIK', time_bpik, resid_bpik);
   %print_stuff('URV - Haar', time_haar, resid_haar);
   %print_stuff('URV - ROS', time_mix, resid_mix);

   times_mean = mean(times_,1);
   times_std = std(times_,[],1);

   print_stuff = @(name,tm,ts,r,n) fprintf('%s\n  time = %1.3e+-%1.3e sec\n  resid = %1.3e\n  norm = %1.3e\n',name,tm,ts,r,n);
   print_stuff('QR', times_mean(1), times_std(1), resid_qr, norm_qr);
   print_stuff('QRCP', times_mean(2), times_std(2), resid_qrcp, norm_qrcp);
   print_stuff('BLENDENPIK', times_mean(3), times_std(3), resid_bpik, norm_bpik);
   print_stuff('URV - Haar', times_mean(4), times_std(4), resid_haar, norm_haar);
   print_stuff('URV - ROS', times_mean(5), times_std(5), resid_mix, norm_mix);

end




function [x] = ls_pgm(A,b)
   function [Q,R] = qr_pgm_block(A)
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

   [m,n] = size(A);
   p = min(m,n);
 
   [Q,R] = qr_pgm_block(A)
   x = backsub(R(1:p,1:p), Q'*b);
   if n > m % underdetermined
      x = [x; zeros(n-m,size(x,2))]; % basic solution, not minimum norm solution
   end

end

function [x] = ls_full_sketch(A,b,p)
   [m,n] = size(A);
   if nargin < 3, p = ceil(.1*m); end
   %TODO info on "counting sketch"

end

function [x] = ls_partial_sketch(A,b,p)
   %TODO
end

function [x] = ls_blendenpik(A,b)
   params = struct;
   params.type='DCT';
   x = blendenpik(A,b,params);
end




function [Q,R] = qr_no_piv(A)
   [m,n] = size(A);
   [Q,R] = qr(A,0);
   % this is another way to measure how far ran(Q) is from ran(A).
   %M = rref([A Q]);
   %nrm = norm(tril(M(:,n+1:end),-1),'fro');
   %ang = subspace(A,Q);
end

function [Q,R,p] = qr_piv(A)
   % this will not work for gpuArray
   [m,n] = size(A);
   [Q,R,p] = qr(A,0);
end

function [Q,R,p] = qr_naive_piv(A)
   % pre-order A based on initial column 2-norms
   [m,n] = size(A);
   [~,p] = sort(sum(A.*conj(A),1), 'descend');
   [Q,R] = qr(A(:,p),0);
end

function [Q,R] = qr_haar(A)
   % sample from the Haar ensemble, which is totally impractical
   % but, we can also sample only once per size of A, which totally biases the timing
   [m,n] = size(A);
   %persistent V
   V = [];
   if isempty(V) || size(V,2) ~= size(A,2)
      [V,~] = qr(randn(size(A,2)),0);
   end
   [Q,R] = qr(A*V',0);
end

