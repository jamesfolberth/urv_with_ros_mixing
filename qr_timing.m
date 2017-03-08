function [] = qr_timing(t_arr_filename, t_arr_filename2)
% QR_TIMING  Timing experiments with QR/QRCP/RURV_ROS
%
% We set up parameters of the timing experiment by modifying the scale_timing function.
% We then run `qr_timing` with no arguments to run the timing experiment.
% To plot a single timing experiment, we do `qr_timing time_file`.
% To plot two timing experiments simultaneously, we do `qr_timing time_file_1 time_file_2`

if nargin == 1
   plot_scaling(t_arr_filename);
   return;
elseif nargin == 2
   plot_scalings(t_arr_filename, t_arr_filename2);
   return;
end

scale_timing();

end

function [] = scale_timing()
% SCALE_TIMING  Run various QR-type factorizations and record runtimes
   now = datetime('now');
   cpu_info = get_cpu_info();

   function [A] = gen_mat(m,n)
      rng(271828);
      %% From Stephen's initial experiments
      %A = randn(m,n) + exp(10*rand(m,n));
      %% If we do this, then more obvious difference///
      %A = bsxfun( @times, A, exp(2*rand(1,n)) );
      %A = qr_gallery('kahan', n);
      A = qr_gallery('condk', [m n], 1e6);
   end
   
   n_sizes = 25;
   n_samp = 5; %10
   n_vec = floor(logspace(2,4,n_sizes)); % 3.5, 4
   m_vec = 2*n_vec; % 10*n_vec

   t_arr_methods = {'DGEQRF', 'DGEQP3', 'RURV_ROS'};
   t_arr = zeros(n_samp, n_sizes, numel(t_arr_methods));

   for i=1:numel(m_vec)
      m = m_vec(i); n = n_vec(i);
      fprintf('n = %d\n', n);
      A = gen_mat(m,n);
      
      for m=1:numel(t_arr_methods)
         [tot, avg, t_arr(:,i,m)] = run_timing(str2func(t_arr_methods{m}), A, n_samp);
      end
   end
   
   t_arr_filename = sprintf('qr_timing_%d%02d%02d_%02d%02d%02d.mat', now.Year, now.Month, now.Day,...
      now.Hour, now.Minute, round(now.Second));
   n_threads = maxNumCompThreads();
   save(t_arr_filename, 't_arr', 't_arr_methods', 'm_vec', 'n_vec', 'cpu_info', 'n_threads');

end

function [] = plot_scaling(t_arr_filename)
% PLOT_SCALING  Plot a single timing experiment
   S = load(t_arr_filename);
   t_arr = S.t_arr;
   t_arr_methods = S.t_arr_methods;
   m_vec = S.m_vec; n_vec = S.n_vec;
   cpu_info = S.cpu_info;
   n_threads = S.n_threads; %JMF: added, will break old runs

   [n_samp, n_sizes, n_methods] = size(t_arr);
   assert(n_methods == numel(t_arr_methods));
   
   means = zeros(n_sizes,n_methods);
   stds = zeros(n_sizes,n_methods);

   for i=1:n_methods
      means(:,i) = mean(t_arr(:,:,i),1);
      stds(:,i) = std(t_arr(:,:,i));
   end
   
   figure(1);
   clf;
   hold on;
   for i=1:n_methods
      loglog(n_vec, means(:,i),'-o', 'linewidth', 2, 'markersize', 2);
      %errorbar(n_vec, means(:,i), stds(:,i)); % dang, this is ugly
   end
   hold off
   axis([min(n_vec) max(n_vec) min(means(:))-max(stds(:)) max(means(:))+max(stds(:))])
   
   if n_threads == 1
      title(sprintf('Average Factorization Runtimes (%d thread)', n_threads), 'fontsize', 20);
   else
      title(sprintf('Average Factorization Runtimes (%d threads)', n_threads), 'fontsize', 20);
   end
   set(gca, 'fontsize', 20);
   %xlabel('m = 2n');
   xlabel('n = m/2');
   ylabel('Average Runtime (s)');
   legend(gca, deal(t_arr_methods), 'Location', 'NorthWest', 'interpreter', 'none');
   set(gca, 'xscale', 'log', 'yscale', 'log');
   axis([1e2 1e4 1e-4 1e3]);
   fprintf(1, cpu_info);

end

function [] = plot_scalings(t_arr_filename, t_arr_filename2)
% PLOT_SCALINGS  Plot two timing experiments on the same set of axes
   S = load(t_arr_filename);
   t_arr = S.t_arr;
   t_arr_methods = S.t_arr_methods;
   m_vec = S.m_vec; n_vec = S.n_vec;
   cpu_info = S.cpu_info;
   n_threads = S.n_threads; %JMF: added, will break old runs

   S2 = load(t_arr_filename2);
   t_arr2 = S2.t_arr;
   t_arr_methods2 = S2.t_arr_methods;
   m_vec2 = S2.m_vec; n_vec2 = S2.n_vec;
   cpu_info2 = S2.cpu_info;
   n_threads2 = S2.n_threads; %JMF: added, will break old runs

   [n_samp, n_sizes, n_methods] = size(t_arr);
   assert(n_methods == numel(t_arr_methods));
   [n_samp2, n_sizes2, n_methods2] = size(t_arr2);
   assert(n_methods2 == numel(t_arr_methods2));
 
   means = zeros(n_sizes,n_methods);
   stds = zeros(n_sizes,n_methods);
   for i=1:n_methods
      means(:,i) = mean(t_arr(:,:,i),1);
      stds(:,i) = std(t_arr(:,:,i));
   end

   means2 = zeros(n_sizes2,n_methods2);
   stds2 = zeros(n_sizes2,n_methods2);
   for i=1:n_methods2
      means2(:,i) = mean(t_arr2(:,:,i),1);
      stds2(:,i) = std(t_arr2(:,:,i));
   end
   
   % Hack to only plot QR,QRCP
   n_methods = n_methods-1;
   n_methods2 = n_methods2-1;
   t_arr_methods = t_arr_methods(1:2);
   t_arr_methods2 = t_arr_methods2(1:2);
 
   figure(1);
   co = get(gca, 'colororder');
   clf;
   hold on;
   for i=1:n_methods
      loglog(n_vec, means(:,i),'-o', 'color', co(i,:), 'linewidth', 2, 'markersize', 10);
      %errorbar(n_vec, means(:,i), stds(:,i), 'color', co(i,:)); % dang, this is ugly
   end
   for i=1:n_methods2
      loglog(n_vec2, means2(:,i), '--x', 'color', co(i,:), 'linewidth', 2, 'markersize', 10);
      %errorbar(n_vec2, means2(:,i), stds2(:,i), 'color', co(i,:)); % dang, this is ugly
   end
   hold off
   %axis([min(n_vec) max(n_vec) min(means(:))-max(stds(:)) max(means(:))+max(stds(:))])
   %axis([min(min(n_vec),min(n_vec2)) max(max(n_vec),max(n_vec2))
   
   for i=1:numel(t_arr_methods)
      if n_threads > 1
         t_arr_methods{i} = sprintf('%s - %d threads', t_arr_methods{i}, n_threads);
      else
         t_arr_methods{i} = sprintf('%s - %d thread', t_arr_methods{i}, n_threads);
      end
   end
   for i=1:numel(t_arr_methods2)
      if n_threads2 > 1
         t_arr_methods2{i} = sprintf('%s - %d threads', t_arr_methods2{i}, n_threads2);
      else
         t_arr_methods2{i} = sprintf('%s - %d thread', t_arr_methods2{i}, n_threads2);
      end
   end
   
   % we don't show error bars, so drop "Average"?
   %title(sprintf('Average Factorization Runtimes', n_threads), 'fontsize', 20);
   title(sprintf('Factorization Runtimes', n_threads), 'fontsize', 20);
   set(gca, 'fontsize', 20);
   %xlabel('m = 2n');
   xlabel('n = m/2');
   ylabel('Average Runtime (s)');
   t_arr_methods
   t_arr_methods2
   legend(gca, deal(horzcat(t_arr_methods,t_arr_methods2)), 'Location', 'NorthWest', 'interpreter', 'none');
   set(gca, 'xscale', 'log', 'yscale', 'log');
   fprintf(1, cpu_info);

end

function [t_tot, t_avg, t_arr] = run_timing(qr_func, A, n_samp, verbosity)
% RUN_TIMING  A driver to run timing experiments for a QR-type functions
   if nargin < 3, n_samp = 10; end
   if nargin < 4, verb = 1; end
      
   % once to JIT compile; hopefully that warms up the CPU governor
   if verb > 0, fprintf('run_timing: Warming up...'); end
   for i=1:1
      qr_func(A);
   end
   if verb > 0, fprintf(' done.\n'); end
   
   
   % now do the actual timing
   t_arr = zeros(n_samp,1);
   t0 = tic();
   for i=1:n_samp
      if verb > 0, fprintf('\rrun_timing: Iteration %d of %d.', i, n_samp); end
      t1 = tic();
      qr_func(A);
      t_arr(i) = toc(t1);
   end
   t_tot = toc(t0);
   t_avg = mean(t_arr);
   if verb > 0, fprintf('\n'); end

end

function [text] = get_cpu_info()
   text = 'CPU info not found.\n';
   if isunix()
      res = system('cat /proc/cpuinfo | head -n 26 > tmpfile.txt');
      if res == 0
         text = fileread('tmpfile.txt');
         delete('tmpfile.txt');
      else
         warning('Not getting CPU info: system call returned with exit status %d', res);
      end
   else
      warning('Not getting CPU info.  This isn''t a UNIX system; I don''t know this.');
   end
end


% Local versions that only do the QR factorization without buliding Q
function [] = DGEQRF(A)
   qr_fact_only(A,0);
end

function [] = DGEQP3(A)
   qr_fact_only(A,1);
end

function [] = RURV_ROS(A)
   
   mix_opt = mix_fftw();
   mix_opt.transpose=1;
   mix_opt.side='r';
   mix_opt.plan_rigor=4;
   mix_opt.wisdom_file='/scratch/wisdom/wisdom.fftw';
   mix_opt.n_threads=maxNumCompThreads(); % JMF: this is deprecated, but it's convenient...
   mix_opt.block_size=250;

   A_mix = mix_fftw(A,1,mix_opt);

   [~,p] = sort(sum(A_mix.*conj(A_mix),1), 'descend'); % post-mixing naive pre-order
   A_mix = A_mix(:,p); 
 
   qr_fact_only(A_mix,0);

end
