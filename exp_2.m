% Filterbank performance on Minnesota traffic graph

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

tic
switch design
    case 'orth'
        fprintf('Designing filters...\n');
        [H0, H1] = orth_design(L, max_eval, 8);
        H = [H0; H1;];        
        toc
        
        fprintf('Determining sampling sets...\n');
        S = set_greedy_deterministic(H, H);
        toc
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(H'*H - 2*eye(n),'fro')^2);
        T = H(S,:)'*H(S,:);
        fprintf('fro_norm_sq(T - I) = %f\n', norm(T - eye(n),'fro')^2);
        toc
    case 'biorth'
        fprintf('Designing filters...\n');
        [H0, H1, G0, G1] = biorth_design(L,max_eval,6,6);
        H = [H0; H1;];
        G = [G0; G1;];
        toc
        
        fprintf('Determining sampling sets...\n');
        S = set_greedy_deterministic(G, H);
        toc
        
        fprintf('fro_norm_sq(G^T*H - I) = %f\n', norm(G'*H - 2*eye(n),'fro')^2);
        T = G(S,:)'*H(S,:);
        fprintf('fro_norm_sq(T - I) = %f\n', norm(T - eye(n),'fro')^2);
        toc
end

T_spec = U'*T*U;

%% plotting

figure;
gplot(A,xy,'k-');
S0 = S(1:n);
S1 = S(n+1:2*n);
xlim([min(xy(:,1))-1 max(xy(:,1))+1]);
ylim([min(xy(:,2))-1 max(xy(:,2))+1]);
% export_fig(['results/graph_1.pdf'],'-transparent');
hold on;
scatter(xy(S0,1),xy(S0,2),20,[1 0 0],'o','Filled');
scatter(xy(S1,1),xy(S1,2),15,[0 0 1],'s','Filled');
axis off;
export_fig(['plots/minn_' design '_sampling_sets.pdf'],'-transparent');

font_size = 20;

figure1 = figure;
scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.35]);
plot(eval, diag(T_spec), 'r. ');
xlim([0 2]);
ylim([0 1.5]);
xlabel('\lambda','FontSize',font_size+4);
ylabel('|T(\lambda)|','FontSize',font_size+4);
set(gca,'FontSize',font_size);
export_fig(['plots/minn_' design '_response.pdf'],'-transparent');

figure2 = figure;
scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure2, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.35]);
plot(eval, max(abs(T_spec - diag(diag(T_spec))), [], 2), 'r. ');
xlim([0 2]);
ylim([0 1.5]);
xlabel('\lambda','FontSize',font_size+4);
ylabel('max_{\mu \neq \lambda} |T(\mu)|','FontSize',font_size+4);
set(gca,'FontSize',font_size);
export_fig(['plots/minn_' design '_leakage.pdf'],'-transparent');
