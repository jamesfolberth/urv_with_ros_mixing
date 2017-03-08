function [] = qlp_test()

if ~add_rrqr_comparison(), return; end % for DGEQPX

%plot_rlvals()
plot_avg_lvals()

end


% In Stewart's QLP approximation paper, the pivoted QLP uses two QRCPs.  What should we be doing?
% JMF: Stewart mentions that in the examples he shows, the 2nd QR doesn't
%      really need to be pivoted; it's mostly for cosmetics (i.e., monotonically 
%      decreasing L-values).  However, for some problems, a 2nd pivoted QR might
%      be necessary.
function [] = plot_rlvals()
   %rng(271828);

   [A,s] = qr_gallery('devil', 128, 22);
   [A,s] = qr_gallery('devil', 128, 22, .9); % all do well for large jumps; smaller are harder
   %A = qr_gallery('kahan', 100); s = svd(A);
   %[A,s] = qr_gallery('rankk', 1000,500);
   %[A,s] = qr_gallery('rankk', 128, 64, 1e-3); % similar to one of Stewart's tests in his QLP paper
   
   cond(A)
   
   %% QR (unpivoted)
   [Q,R] = qr(A,0);
   rval_qr = abs(diag(R));
   rval_rat_qr = abs(diag(R))./s;
   [Q,Lt] = qr(R',0); % Stewart states pivoting is _maybe_ needed
   lval_qr = abs(diag(Lt));

   %% QRCP
   [Q,R,e] = qr(A,'vector');
   rval_qrcp = abs(diag(R));
   rval_rat_qrcp = abs(diag(R))./s;
   %[Q,Lt,e] = qr(R','vector');
   [Q,Lt] = qr(R',0); % Stewart states pivoting is _maybe_ needed
   %[Q,Lt,~] = qr(R',0); 
   lval_qrcp = abs(diag(Lt));

   %% RURV_ROS
   mix_opt = mix_fftw();
   mix_opt.transpose=1;
   mix_opt.side='r';
   [A_mix,darr] = mix_fftw(A,1,mix_opt);
   % naive preorder
   [~,p] = sort(sum(A_mix.*conj(A_mix),1), 'descend'); % post-mixing naive pre-order
	pinv(p) = 1:numel(p);
   A_mix = A_mix(:,p); 
   % unpivoted QR
   [Q,R] = qr(A_mix,0);
   rval_mix = abs(diag(R));
   rval_rat_mix = abs(diag(R))./s;
   [Q,Lt] = qr(R',0); 
   lval_mix = abs(diag(Lt));

   %% Haar mixed QR (very similar to our row mixed QR)
   % with a hybrid RURV_Haar+QRCP QLP, we kinda see some gaps; it's not good, but kinda there
   [V,~] = qr(randn(size(A)),0);
   %[Q,R] = qr(A*V',0);
   
   A_mix = A*V';
   % naive preorder
   [~,p] = sort(sum(A_mix.*conj(A_mix),1), 'descend'); % post-mixing naive pre-order
	pinv(p) = 1:numel(p);
   A_mix = A_mix(:,p); 
   % unpivoted QR
   [Q,R] = qr(A_mix,0);
   
   rval_haar = abs(diag(R));
   rval_rat_haar = abs(diag(R))./s;
   [Q,Lt] = qr(R',0); % unpivoted QR produces better L-values than RURV_Haar
   %[Q,Lt,~] = qr(R',0);
   %[V2,~] = qr(randn(size(A)),0);
   %[U2,Lt] = qr(R'*V2',0);
   lval_haar = abs(diag(Lt));
      
   % More QLP iterations
   % pretty much no benefit
   %for i=1:100
   %   [V,~] = qr(randn(size(Lt')),0);
   %   A_mix = Lt*V';
   %   % naive preorder
   %   [~,p] = sort(sum(A_mix.*conj(A_mix),1), 'descend'); % post-mixing naive pre-order
	%   pinv(p) = 1:numel(p);
   %   A_mix = A_mix(:,p); 
   %   % unpivoted QR
   %   [Q,R] = qr(A_mix,0);
   %   [Q,Lt] = qr(R',0); % unpivoted QR produces better L-values than RURV_Haar
   %   lval_haar = abs(diag(Lt));
   %end
   
   
   %% PGM's blocked QR (is NOT RR)
   % Performs just a bit worse than QRCP+QLP
   [Q,R] = qr_pgm_block(A);
   rval_pgm = abs(diag(R));
   rval_rat_pgm = abs(diag(R))./s;
   [Q,Lt] = qr(R',0);
   %[Q,Lt] = qr_pgm_block(R');
   lval_pgm = abs(diag(Lt));
   
   %% DGEQPX - known RRQR
   % for qr_gallery('devil', 128, 22), this finds the first gap, but then starts to break down
   [~,R,~] = qrxp(A);
   %[Q,R,e] = qrxp(A,0,0,0,size(A,1),2,s);
   rval_dgeqpx = abs(diag(R));
   rval_rat_dgeqpx = abs(diag(R))./s;
   [Q,Lt,~] = qrxp(R',0);
   lval_dgeqpx = abs(diag(Lt));

   %names = {'SVD', 'QRCP', 'QR', 'Row mixed QR', 'PGM Block QR'};
   names = {'SVD', 'QR', 'QRCP', 'RURV\_Haar', 'RURV\_ROS'};

   %figure(1); clf;
   %hold on
   %   plot(s,'k-');
   %   plot(rval_qr,'kd');
   %   plot(rval_qrcp,'go');
   %   plot(rval_haar,'m+:');
   %   plot(rval_mix,'r*:');
   %   %plot(rval_dgeqpx,'bx:');
   %   %plot(rval_pgm,'cs:');
   %hold off
   %title('R-values', 'fontsize', 20);
   %legend(deal(names),'location','Southwest')
   %set(gca,'yscale','log')   
   %xlabel('Column No. i')
   %ylabel('R-values and singular values')
   %set(gca, 'fontsize', 20);

   %figure(2); clf;
   %hold on
   %   plot(rval_rat_qr,'kd');
   %   plot(rval_rat_qrcp,'go');
   %   plot(rval_rat_haar,'m+:');
   %   plot(rval_rat_mix,'r*:');
   %   %plot(rval_rat_dgeqpx,'bx:');
   %   %plot(rval_rat_pgm,'cs:');
   %hold off
   %title('R-values Ratios', 'fontsize', 20);
   %legend(names{2:end},'location','Southwest')
   %set(gca,'yscale','log')   
   %xlabel('Column No. i')
   %ylabel('Ratios |R(i,i)|/\sigma(i)', 'interpreter', 'tex')
   %set(gca, 'fontsize', 20);

   figure(3); clf;
   hold on
      plot(s,'k-');
      plot(lval_qr,'kd');
      plot(lval_qrcp,'go');
      plot(lval_haar,'m+:');
      plot(lval_mix,'r*:');
      %plot(lval_dgeqpx,'bx:');
      %plot(lval_pgm,'cs:');
   hold off
   title('L-values', 'fontsize', 20);
   names_qlp = {names{1}};
   for i=2:numel(names)
      names_qlp{i} = [names{i}, ' + QLP'];
   end
   legend(deal(names_qlp),'location','Southwest')
   set(gca,'yscale','log')   
   xlabel('Column No. i')
   ylabel('L-values and singular values')
   set(gca, 'fontsize', 20);


end


function [] = plot_avg_lvals()
   %rng(271828);

   [A,s] = qr_gallery('devil', 128, 22);
   [A,s] = qr_gallery('devil', 128, 22, .9); % all do well for large jumps; smaller are harder
   %A = qr_gallery('kahan', 100); s = svd(A);
   %[A,s] = qr_gallery('rankk', 1000,500);
   %[A,s] = qr_gallery('rankk', 128, 64, 1e-3); % similar to one of Stewart's tests in his QLP paper
   
   cond(A)
   
   n_samples = 25;
   Lval_qr =     zeros(numel(s), n_samples);
   Lval_qrcp =   zeros(numel(s), n_samples);
   Lval_mix =    zeros(numel(s), n_samples);
   Lval_haar =   zeros(numel(s), n_samples);
   Lval_dgeqpx = zeros(numel(s), n_samples);
   Lval_pgm =    zeros(numel(s), n_samples);
   
   for i=1:n_samples
      %NOTE: We prefer using a different random matrix for each sample
      %      Plot QR+QLP and QRCP+QLP with error bars.
      [A,s] = qr_gallery('devil', 128, 22, .9); % all do well for large jumps; smaller are harder
      %A = qr_gallery('kahan', 128, acos(0.01));
      %s = dgejsv(A);
     

      %% QR (unpivoted)
      [Q,R] = qr(A,0);
      rval_qr = abs(diag(R));
      rval_rat_qr = abs(diag(R))./s;
      [Q,Lt] = qr(R',0); % Stewart states pivoting is _maybe_ needed
      lval_qr = abs(diag(Lt));
      Lval_qr(:,i) = abs(diag(Lt));

      %% QRCP
      [Q,R,e] = qr(A,'vector');
      rval_qrcp = abs(diag(R));
      rval_rat_qrcp = abs(diag(R))./s;
      %[Q,Lt,e] = qr(R','vector');
      [Q,Lt] = qr(R',0); % Stewart states pivoting is _maybe_ needed
      %[Q,Lt,~] = qr(R',0); 
      lval_qrcp = abs(diag(Lt));
      Lval_qrcp(:,i) = abs(diag(Lt));

      %% RURV_ROS
      mix_opt = mix_fftw();
      mix_opt.transpose=1;
      mix_opt.side='r';
      [A_mix,darr] = mix_fftw(A,1,mix_opt);
      % naive preorder
      [~,p] = sort(sum(A_mix.*conj(A_mix),1), 'descend'); % post-mixing naive pre-order
	   pinv(p) = 1:numel(p);
      A_mix = A_mix(:,p); 
      % unpivoted QR
      [Q,R] = qr(A_mix,0);
      rval_mix = abs(diag(R));
      rval_rat_mix = abs(diag(R))./s;
      [Q,Lt] = qr(R',0); 
      lval_mix = abs(diag(Lt));
      Lval_mix(:,i) = abs(diag(Lt));

      %% Haar mixed QR (very similar to our row mixed QR)
      % with a hybrid RURV_Haar+QRCP QLP, we kinda see some gaps; it's not good, but kinda there
      [V,~] = qr(randn(size(A)),0);
      %[Q,R] = qr(A*V',0);
      
      A_mix = A*V';
      % naive preorder
      [~,p] = sort(sum(A_mix.*conj(A_mix),1), 'descend'); % post-mixing naive pre-order
	   pinv(p) = 1:numel(p);
      A_mix = A_mix(:,p); 
      % unpivoted QR
      [Q,R] = qr(A_mix,0);
      
      rval_haar = abs(diag(R));
      rval_rat_haar = abs(diag(R))./s;
      [Q,Lt] = qr(R',0); % unpivoted QR produces better L-values than RURV_Haar
      %[Q,Lt,~] = qr(R',0);
      %[V2,~] = qr(randn(size(A)),0);
      %[U2,Lt] = qr(R'*V2',0);
      lval_haar = abs(diag(Lt));
      Lval_haar(:,i) = abs(diag(Lt));
         
      % More QLP iterations
      % pretty much no benefit
      %for i=1:100
      %   [V,~] = qr(randn(size(Lt')),0);
      %   A_mix = Lt*V';
      %   % naive preorder
      %   [~,p] = sort(sum(A_mix.*conj(A_mix),1), 'descend'); % post-mixing naive pre-order
	   %   pinv(p) = 1:numel(p);
      %   A_mix = A_mix(:,p); 
      %   % unpivoted QR
      %   [Q,R] = qr(A_mix,0);
      %   [Q,Lt] = qr(R',0); % unpivoted QR produces better L-values than RURV_Haar
      %   lval_haar = abs(diag(Lt));
      %end
      
      
      %% PGM's blocked QR (is NOT RR)
      % Performs just a bit worse than QRCP+QLP
      [Q,R] = qr_pgm_block(A);
      rval_pgm = abs(diag(R));
      rval_rat_pgm = abs(diag(R))./s;
      [Q,Lt] = qr(R',0);
      %[Q,Lt] = qr_pgm_block(R');
      lval_pgm = abs(diag(Lt));
      Lval_pgm(:,i) = abs(diag(Lt));
      
      %% DGEQPX - known RRQR
      % for qr_gallery('devil', 128, 22), this finds the first gap, but then starts to break down
      [~,R,~] = qrxp(A);
      %[Q,R,e] = qrxp(A,0,0,0,size(A,1),2,s);
      rval_dgeqpx = abs(diag(R));
      rval_rat_dgeqpx = abs(diag(R))./s;
      [Q,Lt,~] = qrxp(R',0);
      lval_dgeqpx = abs(diag(Lt));
      Lval_dgeqpx(:,i) = abs(diag(Lt));

   end % for i=1:n_samples
   
   %names = {'SVD', 'QR', 'QRCP', 'RURV\_Haar', 'RURV\_ROS'};
   names_mats = {{'SVD', 's'}, {'QR', 'lval_qr'}, {'QRCP', 'lval_qrcp'},...
                 {'RURV\_Haar', 'lval_haar'}, {'RURV\_ROS', 'lval_mix'}};

   lval_qr_mean =     mean(Lval_qr,2);
   lval_qrcp_mean =   mean(Lval_qrcp,2);
   lval_mix_mean =    mean(Lval_mix,2);
   lval_haar_mean =   mean(Lval_haar,2);
   lval_pgm_mean =    mean(Lval_pgm,2);
   lval_dgeqpx_mean = mean(Lval_dgeqpx,2);

   lval_qr_median =     median(Lval_qr,2);
   lval_qrcp_median =   median(Lval_qrcp,2);
   lval_mix_median =    median(Lval_mix,2);
   lval_haar_median =   median(Lval_haar,2);
   lval_pgm_median =    median(Lval_pgm,2);
   lval_dgeqpx_median = median(Lval_dgeqpx,2);

   lval_qr_std =     std(Lval_qr,[],2);
   lval_qrcp_std =   std(Lval_qrcp,[],2);
   lval_mix_std =    std(Lval_mix,[],2);
   lval_haar_std =   std(Lval_haar,[],2);
   lval_pgm_std =    std(Lval_pgm,[],2);
   lval_dgeqpx_std = std(Lval_dgeqpx,[],2);

   lval_qr_min =     min(Lval_qr,[],2);
   lval_qrcp_min =   min(Lval_qrcp,[],2);
   lval_mix_min =    min(Lval_mix,[],2);
   lval_haar_min =   min(Lval_haar,[],2);
   lval_pgm_min =    min(Lval_pgm,[],2);
   lval_dgeqpx_min = min(Lval_dgeqpx,[],2);
   min_all = min([lval_qr_min; lval_qrcp_min; lval_mix_min; lval_haar_min]);
 
   lval_qr_max =     max(Lval_qr,[],2);
   lval_qrcp_max =   max(Lval_qrcp,[],2);
   lval_mix_max =    max(Lval_mix,[],2);
   lval_haar_max =   max(Lval_haar,[],2);
   lval_pgm_max =    max(Lval_pgm,[],2);
   lval_dgeqpx_max = max(Lval_dgeqpx,[],2);
   max_all = max([lval_qr_max; lval_qrcp_max; lval_mix_max; lval_haar_max]);
 
   holdMarker('reset');
   for i=2:numel(names_mats)
      figure(i-1); clf;
      marker = holdMarker();
      name = names_mats{i}{1};
      basename = names_mats{i}{2};
      means = eval(strcat(basename, '_mean'));
      medians = eval(strcat(basename, '_median'));
      mins = eval(strcat(basename, '_min'));
      maxs = eval(strcat(basename, '_max'));
      stds = eval(strcat(basename,'_std'));
      hold on
      plot(s, 'k-');
      %if strcmp(name, 'QR') || strcmp(name, 'QRCP') % not random, so no errorbar
      %   plot(means, 'bd');
      %else
         %errorbar(means, stds, 'bd'); % can't change the width, and all hacks I tried failed
			plot(means, 'bd');
			myerrorbar(1:numel(s), means, mins, maxs, 0.5, 'b');
      %end
      hold off
      title(strcat({'L-values - '}, {name}), 'fontsize', 20);
      ylabel('L-values');
      xlabel('Column No. i');
      set(gca, 'fontsize', 20);
      set(gca, 'yscale', 'log');
      %axis([1 128 1e-5 1e0]);
      axis([1 128 min(1e-5, min_all) max(max_all,1e0)]);
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


