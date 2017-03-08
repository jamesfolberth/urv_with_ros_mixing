function [] = column_norm_smoothing()

rng(271828)

m = 250; n = 250;
sort_norms = 0;
A = randn(m,n);

% From Stephen's initial experiments
A = randn(m,n) + exp(10*rand(m,n));
% If we do this, then more obvious difference
A = bsxfun( @times, A, exp(2*rand(1,n)) );
A = A/mean(norms(A));



% sort by column norms
nrms = norms(A);
if sort_norms, [nrms,p] = sort(nrms,2,'ascend'); end


[V,~] = qr(randn(n));
nrms_haar = norms(A*V');
if sort_norms, [nrms_haar,p] = sort(nrms_haar,2,'ascend'); end


mix_opt = mix_fftw();
mix_opt.transpose=1;
mix_opt.side='r';
A_mix = mix_fftw(A,1,mix_opt);
nrms_mix = norms(A_mix);
if sort_norms, [nrms_mix,p] = sort(nrms_mix,2,'ascend'); end


figure(1); clf;
%suptitle('Effect of Mixing on Column Norms of A', 'fontsize', 20);
subplot(1,2,1);
set(gca, 'fontsize', 20);
set(gcf, 'position', [2 23 1918 800]); % make it a bit shorter+fatter
co = get(gca, 'colororder');
hold on;
plot(1:n,nrms,'o','Color', co(1,:), 'MarkerFaceColor', co(1,:));
plot(1:n,nrms_haar,'o', 'Color', co(4,:), 'MarkerFaceColor', co(4,:));
plot(1:n,nrms_mix,'o','Color', co(3,:), 'MarkerFaceColor', co(3,:));
hold off;
%set(gca, 'yscale', 'log')

%if sort_norms
%   title('Effect of Mixing on Column Norms of A - Sorted');
%   xlabel('Sorted Column Index');
%else
%   title('Effect of Mixing on Column Norms of A');
%   xlabel('Column Index');
%end
legend('A - No mixing', 'A - Haar mixing', 'A - ROS mixing', 'Location', 'NorthWest');
ylabel('Column 2-Norm');
xlabel('Column Index');
axis([1 n 0 2.5])


subplot(1,2,2);
set(gca, 'fontsize', 20);
hold on;
boxplot([nrms.' nrms_haar.' nrms_mix.'], {'No mixing', 'Haar mixing', 'ROS mixing'});
hold off;
%title('Effect of Mixing on Column Norms of A')
%ylabel('Column 2-Norm');
v = axis;
axis([v(1) v(2) 0 2.5])

% set box colors
% https://www.mathworks.com/matlabcentral/answers/122426-change-color-of-a-box-in-boxplot
a = get(get(gca, 'children'), 'children');
set(a(7), 'Color', co(3,:)) % in reverse order...
set(a(8), 'Color', co(4,:))
set(a(9), 'Color', co(1,:))

[mean(nrms) std(nrms)]
[mean(nrms_haar) std(nrms_haar)]
[mean(nrms_mix) std(nrms_mix)]

end
