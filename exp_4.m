% Reconstruction error for random signals - biorth design

addpath(genpath('exportfig'));

design = 'biorth';

rng('default');

% % Generate graph
graph = 1;
switch graph
    case 1 % --- Ring graph + random cross-links
        N = 10;
        p = 0.4;
        A = diag(ones(1,N-1),1);
        A(1,end) = 1;
        for i = 1:N-1
            for j = i+1:N
                if A(i,j) == 0
                    A(i,j) = rand(1) < p;
                end
            end
        end
        A = A + A';
        theta = 2*pi/N;
        xy = [cos((0:N-1)'*theta) sin((0:N-1)'*theta)];
    case 2 % --- Erdos-Renyi random graph
        N = 10;
        p = 0.5;
        A = rand(N) < p;
        A = A - tril(A);
        A = A + A';
        xy = rand(N);
end

% % Create Laplacian and compute spectrum

% [num_conn_comp, conn_ind] = graphconncomp(sparse(A));
% comps = unique(conn_ind);
% A = A(conn_ind==comps(1), conn_ind==comps(1));
% xy = xy(conn_ind==comps(1),:);

N = size(A,1);
Deg = diag(sum(A,2));
% L = Deg - A;
L = eye(N) - Deg^(-1/2) * A * Deg^(-1/2);

[U, temp] = eig(L);
[eval,perm] = sort(real(diag(temp)));
U = U(:,perm);

% max_eval = max(eval);
max_eval = 2;

% % create filters
switch design
    case 'orth'
        fprintf('Designing filters...\n');
        [H0, H1] = orth_design(L, max_eval, 8);
        H = [H0; H1;];        
        
        fprintf('Determining sampling sets...\n');
        S_greedy = set_greedy_deterministic(H, H);
        S_opt = set_exhaustive(H, H);
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(H'*H - 2*eye(N),'fro')^2);
        T_greedy = H(S_greedy,:)'*H(S_greedy,:);
        fprintf('fro_norm_sq(T_greedy - I) = %f\n', norm(T_greedy - eye(N),'fro')^2);
        T_opt = H(S_opt,:)'*H(S_opt,:);
        fprintf('fro_norm_sq(T_opt - I) = %f\n', norm(T_opt - eye(N),'fro')^2);
        
    case 'biorth'
        fprintf('Designing filters...\n');
        [H0, H1, G0, G1] = biorth_design(L,max_eval,6,6);
        H = [H0; H1;];
        G = [G0; G1;];
        
        fprintf('Determining sampling sets...\n');
        S_greedy = set_greedy_deterministic(G, H);
        S_opt = set_exhaustive(G, H);
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(G'*H - 2*eye(N),'fro')^2);
        T_greedy = G(S_greedy,:)'*H(S_greedy,:);
        fprintf('fro_norm_sq(T_greedy - I) = %f\n', norm(T_greedy - eye(N),'fro')^2);
        T_opt = G(S_opt,:)'*H(S_opt,:);
        fprintf('fro_norm_sq(T_opt - I) = %f\n', norm(T_opt - eye(N),'fro')^2);
end

%%

num_reps = 100;
num_rand_samp_patterns = 100;
rel_err_random = zeros(num_reps,num_rand_samp_patterns);
rel_err_greedy = zeros(num_reps,1);
rel_err_opt = zeros(num_reps,1);

for i = 1:num_reps
    f = randn(N,1);
    
    for j = 1:num_rand_samp_patterns
        S0 = randsample(1:N,N/2);
        temp = false(N,1);
        temp(S0) = true;
        S_random = [temp; ~temp];
        
        f_tx = H(S_random,:) * f;
        switch design
            case 'orth'
                f_recon = H(S_random,:)' * f_tx;
            case 'biorth'
                f_recon = G(S_random,:)' * f_tx;
        end
        err = f - f_recon;
        rel_err_random(i,j) = (norm(err) / norm(f))^2;
    end
    
    f_tx = H(S_greedy,:) * f;
    switch design
        case 'orth'
            f_recon = H(S_greedy,:)' * f_tx;
        case 'biorth'
            f_recon = G(S_greedy,:)' * f_tx;
    end
    err = f - f_recon;
    rel_err_greedy(i) = (norm(err) / norm(f))^2;
    
    f_tx = H(S_opt,:) * f;
    switch design
        case 'orth'
            f_recon = H(S_opt,:)' * f_tx;
        case 'biorth'
            f_recon = G(S_opt,:)' * f_tx;
    end
    err = f - f_recon;
    rel_err_opt(i) = (norm(err) / norm(f))^2;
end

mean_rel_err_random = mean(rel_err_random(:));
mean_rel_err_greedy = mean(rel_err_greedy);
mean_rel_err_opt = mean(rel_err_opt);

std_rel_err_random = std(rel_err_random(:));
std_rel_err_greedy = std(rel_err_greedy);
std_rel_err_opt = std(rel_err_opt);

fprintf('rel_err_random = %f +- %f\n', mean_rel_err_random, std_rel_err_random);
fprintf('rel_err_greedy = %f +- %f\n', mean_rel_err_greedy, std_rel_err_greedy);
fprintf('rel_err_opt    = %f +- %f\n', mean_rel_err_opt, std_rel_err_opt);

%% plotting

figure;
gplot(A,xy,'k-');
min_x = min(xy(:,1)); max_x = max(xy(:,1)); range_x = max_x - min_x;
min_y = min(xy(:,2)); max_y = max(xy(:,2)); range_y = max_y - min_y;
xlim([min_x-0.1*range_x max_x+0.1*range_x]);
ylim([min_y-0.1*range_y max_y+0.1*range_y]);
hold on;
scatter(xy(:,1),xy(:,2),150,[0 0 0],'o','Filled');
axis equal;
axis off;
export_fig(['plots/graph_ring_p' num2str(100*p) '.pdf'],'-transparent');

figure;
gplot(A,xy,'k-');
S0 = S_greedy(1:N);
S1 = S_greedy(N+1:2*N);
min_x = min(xy(:,1)); max_x = max(xy(:,1)); range_x = max_x - min_x;
min_y = min(xy(:,2)); max_y = max(xy(:,2)); range_y = max_y - min_y;
xlim([min_x-0.1*range_x max_x+0.1*range_x]);
ylim([min_y-0.1*range_y max_y+0.1*range_y]);
hold on;
scatter(xy(~(S0|S1),1),xy(~(S0|S1),2),200,[0 0 0],'o');
scatter(xy(S0,1),xy(S0,2),200,[1 0 0],'v','Filled');
scatter(xy(S1,1),xy(S1,2),200,[0 0 1],'^','Filled');
axis equal;
axis off;
export_fig(['plots/graph_ring_p' num2str(100*p) '_' design '_greedy_samp.pdf'],'-transparent');

figure;
gplot(A,xy,'k-');
S0 = S_opt(1:N);
S1 = S_opt(N+1:2*N);
min_x = min(xy(:,1)); max_x = max(xy(:,1)); range_x = max_x - min_x;
min_y = min(xy(:,2)); max_y = max(xy(:,2)); range_y = max_y - min_y;
xlim([min_x-0.1*range_x max_x+0.1*range_x]);
ylim([min_y-0.1*range_y max_y+0.1*range_y]);
hold on;
scatter(xy(~(S0|S1),1),xy(~(S0|S1),2),200,[0 0 0],'o');
scatter(xy(S0,1),xy(S0,2),200,[1 0 0],'v','Filled');
scatter(xy(S1,1),xy(S1,2),200,[0 0 1],'^','Filled');
axis equal;
axis off;
export_fig(['plots/graph_ring_p' num2str(100*p) '_' design '_opt_samp.pdf'],'-transparent');

