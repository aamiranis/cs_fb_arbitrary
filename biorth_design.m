function [H0, H1, G0, G1] = biorth_design(L, max_eval, k0, k1)

addpath(genpath('graphBior_Filterbanks'));
[h0, h1] = biorth_kernel(k0,k1);

% h0 = h0/sqrt(2);
% h1 = h1/sqrt(2);

scale0 = h0(1);
scale1 = h1(1);

rh0 = roots(h0) * max_eval/2;
rh1 = roots(h1) * max_eval/2;

rg0 = (2 - roots(h1)) * max_eval/2;
rg1 = (2 - roots(h0)) * max_eval/2;

h0 = scale0 * poly(rh0) * (2/max_eval)^length(rh0);
h1 = scale1 * poly(rh1) * (2/max_eval)^length(rh1);

g0 = (-1)^(length(h1)-1) * scale1 * poly(rg0) * (2/max_eval)^length(rg0);
g1 = (-1)^(length(h0)-1) * scale0 * poly(rg1) * (2/max_eval)^length(rg1);

% length(h0)
% length(h1)
% length(g0)
% length(g1)

% t = (0:0.1:max_eval)';
% fh0 = polyval(h0,t);
% fh1 = polyval(h1,t);
% fg0 = polyval(g0,t);
% fg1 = polyval(g1,t);
% 
% plot(t,[fh0 fh1 fg0 fg1 fg0.*fh0 + fg1.*fh1],'LineWidth',2);
% legend('h_0(\lambda)','h_1(\lambda)','g_0(\lambda)','g_1(\lambda)','h_0(\lambda) g_0(\lambda) + h1(\lambda) g_1(\lambda)');
% ylim([-0.5 2.5]);
% xlim([0 max_eval]);
% grid on;
% % export_fig('results/biorth_design.pdf','-transparent');

H0 = polyvalm(h0,L);
H1 = polyvalm(h1,L);

G0 = polyvalm(g0,L);
G1 = polyvalm(g1,L);