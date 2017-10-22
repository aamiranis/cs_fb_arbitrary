function [H0, H1] = orth_design(L, max_eval, filter_length)

addpath(genpath('graphBior_Filterbanks'));

kernel = @(x) meyer_kernel(2*x/max_eval);

freq_range = [0 max_eval];

c0 = sgwt_cheby_coeff(kernel, filter_length, filter_length+1, freq_range);

% t = (0:0.1:max_eval)';
% g0 = sgwt_cheby_eval(t,c0,freq_range);
% g1 = sgwt_cheby_eval(max_eval - t,c0,freq_range);
% figure, plot(t,[g0 g1 g0.*g0 + g1.*g1],'LineWidth',2);
% legend('h_0(\lambda)', 'h_1(\lambda)', 'h_0^2(\lambda) + h_1^2(\lambda)')
% % hold on;
% % plot(t,kernel(t),'k-');
% ylim([-0.5 2.5]);
% xlim([0 max_eval]);
% grid on;
% % export_fig('results/orth_design.pdf','-transparent');


n = length(L);

H0 = zeros(n);
H1 = zeros(n);

for i = 1:length(L)
    e = zeros(n,1);
    e(i) = 1;
    
    H0(:,i) = sgwt_cheby_op(e,L,c0,freq_range);
    H1(:,i) = sgwt_cheby_op(e,max_eval*eye(n) - L, c0,freq_range);
end