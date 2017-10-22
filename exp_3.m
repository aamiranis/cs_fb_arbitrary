% Reconstruction error for random signals

addpath(genpath('exportfig'));

design = 'orth';

% % Generate graph
load('graphBior_Filterbanks/Datasets/min_traffic_graph.mat');
load('graphBior_Filterbanks/Datasets/min_graph_signal.mat');

% % Create Laplacian and compute spectrum

[num_conn_comp, conn_ind] = graphconncomp(sparse(A));
comps = unique(conn_ind);
A = A(conn_ind==comps(1), conn_ind==comps(1));
xy = xy(conn_ind==comps(1),:);

n = size(A,1);
Deg = diag(sum(A,2));
% L = Deg - A;
L = eye(n) - Deg^(-1/2) * A * Deg^(-1/2);

[U, temp] = eig(L);
[eval,perm] = sort(real(diag(temp)));
U = U(:,perm);

% max_eval = max(eval);
max_eval = 2;

% % create filters

% random samples
rng('default');
S0 = randsample(1:n,n/2);
S_random = false(n,1);
S_random(S0) = true;
S_random = [S_random; ~S_random];

% maxcut sampling
S_maxcut = U(:,end) > 0;
S_maxcut = [S_maxcut; ~S_maxcut];

tic
switch design
    case 'orth'
        fprintf('Designing filters...\n');
        [H0, H1] = orth_design(L, max_eval, 8);
        H = [H0; H1;];        
        toc
        
        fprintf('Determining sampling sets...\n');
        S_greedy = set_greedy_deterministic(H, H);
        toc
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(H'*H - 2*eye(n),'fro')^2);
        T_greedy = H(S_greedy,:)'*H(S_greedy,:);
        fprintf('fro_norm_sq(T_greedy - I) = %f\n', norm(T_greedy - eye(n),'fro')^2);
        toc
        
        T_random = H(S_random,:)'*H(S_random,:);
        fprintf('fro_norm_sq(T_random - I) = %f\n', norm(T_random - eye(n),'fro')^2);
        
        T_maxcut = H(S_maxcut,:)'*H(S_maxcut,:);
        fprintf('fro_norm_sq(T_maxcut - I) = %f\n', norm(T_maxcut - eye(n),'fro')^2);
    case 'biorth'
        fprintf('Designing filters...\n');
        [H0, H1, G0, G1] = biorth_design(L,max_eval,6,6);
        H = [H0; H1;];
        G = [G0; G1;];
        toc
        
        fprintf('Determining sampling sets...\n');
        S_greedy = set_greedy_deterministic(G, H);
        toc
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(G'*H - 2*eye(n),'fro')^2);
        T_greedy = G(S_greedy,:)'*H(S_greedy,:);
        fprintf('fro_norm_sq(T_greedy - I) = %f\n', norm(T_greedy - eye(n),'fro')^2);
        toc
        
        T_random = G(S_random,:)'*H(S_random,:);
        fprintf('fro_norm_sq(T_random - I) = %f\n', norm(T_random - eye(n),'fro')^2);
        
        T_maxcut = G(S_maxcut,:)'*H(S_maxcut,:);
        fprintf('fro_norm_sq(T_maxcut - I) = %f\n', norm(T_maxcut - eye(n),'fro')^2);
end

%%

num_reps = 1000;
rel_err_greedy = zeros(num_reps,1);
rel_err_random = zeros(num_reps,1);
rel_err_maxcut = zeros(num_reps,1);

for i = 1:num_reps
    f = randn(n,1);
    
    f_tx = H(S_greedy,:) * f;
    switch design
        case 'orth'
            f_recon = H(S_greedy,:)' * f_tx;
        case 'biorth'
            f_recon = G(S_greedy,:)' * f_tx;
    end
    err = f - f_recon;
    rel_err_greedy(i) = (norm(err) / norm(f))^2;
    
    f_tx = H(S_random,:) * f;
    switch design
        case 'orth'
            f_recon = H(S_random,:)' * f_tx;
        case 'biorth'
            f_recon = G(S_random,:)' * f_tx;
    end
    err = f - f_recon;
    rel_err_random(i) = (norm(err) / norm(f))^2;
    
    f_tx = H(S_maxcut,:) * f;
    switch design
        case 'orth'
            f_recon = H(S_maxcut,:)' * f_tx;
        case 'biorth'
            f_recon = G(S_maxcut,:)' * f_tx;
    end
    err = f - f_recon;
    rel_err_maxcut(i) = (norm(err) / norm(f))^2;
end

mean_rel_err_greedy = mean(rel_err_greedy);
mean_rel_err_random = mean(rel_err_random);
mean_rel_err_maxcut = mean(rel_err_maxcut);

toc

% save('~/tmp/exp_3a.mat');