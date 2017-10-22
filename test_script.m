% addpath(genpath('~/matlab/plotting'));

design = 'biorth';

% % Generate graph

rng('default');

n = 1000;
A = rand(n);
A  = A < 0.01;

A = (A+A')>0;

% % Create Laplacian and compute spectrum

Deg = diag(sum(A,2));
L = Deg - A;
% L = eye(n) - Deg^(-1/2) * A * Deg^(-1/2);

[U, temp] = eig(L);
[eval,perm] = sort(real(diag(temp)));
U = U(:,perm);

max_eval = max(eval);

n = size(A,1);

% % create filters

switch design
    case 'orth'
        [H0, H1] = orth_design(L,max_eval,8);
        H = [H0; H1;];        
        
%         S = set_greedy_submodular(H, H);
%         S = set_spectral(H, H);
%         S = set_greedy(H, H);
        S = set_greedy_randomized(H, H);
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(H'*H - 2*eye(n),'fro')^2);
        T = H(S,:)'*H(S,:);
        fprintf('fro_norm_sq(T - I) = %f\n', norm(T - eye(n),'fro')^2);

    case 'biorth'
        [H0, H1, G0, G1] = biorth_design(L,max_eval,4,4);
        H = [H0; H1;];
        G = [G0; G1;];
        
%         S = set_greedy_submodular(G, H);
%         S = set_spectral(G, H);
%         S = set_greedy(G, H);
        S = set_greedy_randomized(G, H);
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(G'*H - 2*eye(n),'fro')^2);
        T = G(S,:)'*H(S,:);
        fprintf('fro_norm_sq(T - I) = %f\n', norm(T - eye(n),'fro')^2);
end

figure;
plot(eval, diag(U'*T*U), 'k. ');
xlim([0 max_eval]);
ylim([0 1.5]);
xlabel('\lambda');
ylabel('Response (approx.)');
% export_fig(['results/graph_2_' design '_response.pdf'],'-transparent');
