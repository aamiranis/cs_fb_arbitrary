% Testing sampling pattern obtained for grid graph

addpath(genpath('exportfig'));

design = 'biorth';

% % Generate graph


A = [
     0     1     0     0     0     0     0     0     0     1
     1     0     1     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     1
     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     1
     1     0     1     0     0     0     1     0     1     0
    ];

N = size(A,1);

theta = 2*pi/N;
Coordinates = [cos((0:N-1)'*theta) sin((0:N-1)'*theta)];

% % Create Laplacian and compute spectrum

n = size(A,1);
Deg = diag(sum(A,2));
% L = Deg - A;
L = eye(n) - Deg^(-1/2) * A * Deg^(-1/2);

[U, temp] = eig(L);
[eval,perm] = sort(real(diag(temp)));
U = U(:,perm);

max_eval = max(eval);

% % create filters

switch design
    case 'orth'
        [H0, H1] = orth_design(L,2,4);
        H = [H0; H1;];        
        
        S = set_greedy_deterministic(H, H);
%         S = set_greedy_deterministic_disjoint(H, H);
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(H'*H - 2*eye(n),'fro')^2);
        T = H(S,:)'*H(S,:);
        fprintf('fro_norm_sq(T - I) = %f\n', norm(T - eye(n),'fro')^2);

    case 'biorth'
        [H0, H1, G0, G1] = biorth_design(L,max_eval,6,6);
        H = [H0; H1;];
        G = [G0; G1;];
        
        S = set_greedy_deterministic(G, H);
%         S = set_greedy_deterministic_disjoint(G, H);
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(G'*H - 2*eye(n),'fro')^2);
        T = G(S,:)'*H(S,:);
        fprintf('fro_norm_sq(T - I) = %f\n', norm(T - eye(n),'fro')^2);
end


figure;
gplot(A,Coordinates,'k-');
S0 = S(1:n);
S1 = S(n+1:2*n);
set(gca,'LineWidth',3);
xlim([min(Coordinates(:,1))-1 max(Coordinates(:,1))+1]);
ylim([min(Coordinates(:,2))-1 max(Coordinates(:,2))+1]);
% export_fig(['results/graph_1.pdf'],'-transparent');
hold on;
scatter(Coordinates(S0,1),Coordinates(S0,2),100,[1 0 0],'^','Filled');
scatter(Coordinates(S1,1),Coordinates(S1,2),100,[0 0 1],'v','Filled');
scatter(Coordinates(~(S0|S1),1),Coordinates(~(S0|S1),2),100,[0 0 0],'o');
axis equal;
axis off;
export_fig(['plots/graph_3_' design '_sampling_sets.pdf'],'-transparent');
% export_fig(['plots/graph_3_' design '_sampling_sets_disjoint.pdf'],'-transparent');

% figure;
% plot(eval, diag(U'*T*U));
% xlim([0 max_eval]);
% ylim([0 1.5]);
% xlabel('\lambda');
% ylabel('Response (approx.)');
% % export_fig(['results/graph_1_' design '_response.pdf'],'-transparent');
