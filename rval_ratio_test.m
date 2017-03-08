function [] = rval_ratio_test()

if ~add_lapack(), return; end % for interface to LAPACK (DGEJSV)

test1()

end

%TODO matrix #15 "spikes" gives us some trouble; mixing twice helps
function [] = test1()
   rng(271828);   

   n = 256;
   verb = 1;
   
   rmethods = {'QRCP', 'RURV_Haar', 'RURV_ROS - n_mix=1', 'RURV_ROS - n_mix=2'}; 
   rnames = {'QRCP', {'RURV\_Haar';'n\_mix = 1'}, {'RURV\_ROS';'n\_mix = 1'}, {'RURV\_ROS';'n\_mix = 2'}};

   rstats = zeros(18,3,numel(rmethods));
   bounds = zeros(18,2,numel(rmethods));

   for id=1:18
      if verb, fprintf(1,'\rid = %d', id); end
      A = qr_gallery_matrix_by_num(id,n);
      %s = svd(A);
      
      %TODO lower/upper bounds with QRPGM?
      %[Q,R,e] = qr(A,0);
      %A0 = A(:,e);
      %bounds(id,1) = eps(1)*min(norms(A0));
      %bounds(id,2) = eps(1)*max(norms(A0));
      %bounds(id,1)
      %bounds(id,2)

      for mid=1:numel(rmethods)
         if strcmp(rmethods{mid}, 'RURV_ROS - n_mix=1')
            R = RURV_ROS(A,1);
         elseif strcmp(rmethods{mid}, 'RURV_ROS - n_mix=2')
            R = RURV_ROS(A,2);
         elseif strcmp(rmethods{mid}, 'RURV_Haar')
            R = RURV_Haar(A);
         else
            method = str2func(rmethods{mid});
            R = method(A);
         end
         %NOTE: for RURV*, the QRCP bounds are only approximate bounds
         %      the true bounds use the sorted R values, but it's 
         %      interesting to see how we shape up to the QRCP bounds.
         r = diag(R);
         %s = svd(R);
         s = dgejsv(R);
         Yt = bsxfun(@times, R, 1./r);
         bounds(id,1,mid) = 1./norm(Yt);
         %bounds(id,2,mid) = norm(inv(Yt)); % <= norm(Yt)
         bounds(id,2,mid) = max(norm(inv(Yt)), norm(Yt)); % <= norm(Yt)

         rstats(id,:,mid) = [min(abs(r./s)) median(abs(r./s)) max(abs(r./s))];
      end
   end
   if verb, fprintf(1,'\r                          \r'); end

   figure(1);
   clf;
   for ind=1:4
      subplot(4,1,ind);
      hold on;
      bar(squeeze(rstats(:,:,ind)), 'BaseValue', 1);
      ylabel(rnames{ind}) 
      plot_bounds(squeeze(bounds(:,:,ind)));
      common_plot_setup();
      hold off;

      if ind == 3
         common_plot_setup(8);
      end
   end

end

function [] = common_plot_setup(e)
   if nargin < 1, e=4; end

   set(gca, 'yscale', 'log')
   set(gca, 'fontsize', 20)
   axis([0 23 10^-e 10^e])
   h_legend = legend('min', 'median', 'max', 'QRCP lower bnd', 'QRCP upper bnd', 'Location', 'East');
   set(h_legend, 'fontsize', 16);
   set(gca, 'XTick', 1:18);
   set(gca, 'YTick', logspace(-e,e,e+1))
end

function [] = plot_bounds(bounds)
   for i=1:size(bounds,1)
      lb = bounds(i,1); ub = bounds(i,2);
      plot([i-0.5, i+0.5], [lb lb], '--r')
      plot([i-0.5, i+0.5], [ub ub], '--b')
   end
end

function [R] = QRCP(A)
   [Q,R,e] = qr(A,0);
end

function [R] = RURV_Haar(A)
   [V,~] = qr(randn(size(A)),0);
   X = qr(A*V',0);
   R = triu(X);
end

function [R] = RURV_ROS(A,n_mix)
   if nargin < 2, n_mix = 1; end

   mix_opt = mix_fftw();
   mix_opt.transpose=1;
   mix_opt.side='r';
   mix_opt.plan_rigor=4;
   mix_opt.n_threads=1;
   mix_opt.block_size=256;

   A_mix = mix_fftw(A,n_mix,mix_opt);

   [~,p] = sort(sum(A_mix.*conj(A_mix),1), 'descend'); % post-mixing naive pre-order
   A_mix = A_mix(:,p); 
   X = qr(A_mix,0);
   R = triu(X);
end
