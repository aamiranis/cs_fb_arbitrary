function [beta bptG beta_dist Colorednodes] = opt_decomp5(A)

% finds optimal bpt decomp of arbitrary graphs
A1 = A;
N = length(A);

% reweight the graph

wA = 0*A;
edges = sum(sum(A))/2;
% figure,plot(xy(:,1),xy(:,2),'*')
% hold on
theta = 1;
while(edges)
    edges = 0;
    B = A + speye(N);
    [p , ~, r s] = dmperm(B);
    C = A(p,p);
    for i = 1:(length(r)-1)
        Atemp = C(r(i):r(i+1) -1,r(i):r(i+1) -1);
        %         Ntemp = size(Atemp,2);
        d = sum(Atemp,2);
        d_inv = d.^(-0.5);
        if sum(d)
            L = Atemp;
            L = diag(d_inv)*L*diag(d_inv);
            L = 0.5*(L+L');
            %             if theta < 6
            opt.disp = 0 ;                      % turn off printing in eigs
            opt.tol = sqrt (eps) ;
            [U,Lam] = eigs (L, 2, 'SA', opt) ;    % find the Fiedler vector v
            Utemp = U(:,2) >= 0;
            %             else
            %                 [U,Lam] = eigs (L, 1,'LM',opt) ;    % find the max_cut
            %                 Utemp = U >= 0;
            %             end
            
            S1 = r(i) + find(Utemp) - 1;
            S2 = r(i) + find(~Utemp) - 1;
            size1 = length(S1);
            size2 = length(S2);
            factor = 2*size1*size2/(N*sum(sum(A(p(S1),p(S2)))));
            edges = edges + sum(sum(A(p(S1),p(S2))));
            AA = 0*A;
            AA(p(S1),p(S2)) = double(A(p(S1),p(S2))>0)*factor;
            AA(p(S2),p(S1)) = double(A(p(S2),p(S1))>0)*factor;
            wA = wA + AA;
            %                         [X Y] = gplot(AA,xy);
            %                         plot(X,Y,'r-','linewidth',2*factor);
            A(p(S1),p(S2)) = 0;
            A(p(S2),p(S1)) = 0;
        end
    end
    %         pause
    if ~edges
        break;
    end
    theta = theta+1;
end


% run max-cut on wA
theta =1;
Adj = wA;
edges = sum(sum(Adj))/2;
beta= [];
NN = length(Adj);
while(edges)
    edges = 0;
    beta(1:NN,theta) = -1;
    B = Adj + speye(NN);
    [p q r s] = dmperm(B);
    C = Adj(p,p);
    for i = 1:(length(r)-1)
        Atemp = C(r(i):r(i+1) -1,s(i):s(i+1) -1);
        %         Ntemp = size(Atemp,2);
        d = sum(Atemp,2);
        d_inv = d.^(-0.5);
        if sum(d)
            L = diag(d) - Atemp;
            L = diag(d_inv)*L*diag(d_inv);
            L = 0.5*(L+L');
            [U Lam] = eigs(L,1);
            U = U >= 0;
            S1 = r(i) + find(U) - 1;
            S2 = r(i) + find(~U) - 1;
            beta(p(S1),theta) = 1 ;
            beta(p(S2),theta) = -1 ;
            edges = edges + sum(sum(Adj(p(S1),p(S2))));
            Adj(p(S1),p(S2)) = 0;
            Adj(p(S2),p(S1)) = 0;
        end
    end
    if ~edges
        theta = theta -1;
        beta = beta(:,1:theta);
        break;
    end
    theta = theta+1;
end
theta =2;
bptG = zeros(N,N,theta);
% A1 = A;
for i = 1:theta
    S1 = find(beta(:,i) == 1);
    S2 = find(beta(:,i) == -1);
    bptG(S1,S2,i) = A1(S1,S2);
    bptG(S2,S1,i) = A1(S2,S1);
    beta(S1,i) = 1;
    A1 = A1 - bptG(:,:,i);
    %     num_edges(i) = sum(sum(bptG(:,:,i)));
end

edges = sum(sum(bptG,1),2)/2;
edges = edges(:);
cum_edges = cumsum(edges)/sum(edges);
theta1 = find(cum_edges > 0.97, 1,'first');
% theta = theta1;
theta = 2;
bptG = bptG(:,:,1:theta);
beta = beta(:,1:theta);

Colorednodes{1} = find(beta(:,1) ==1 & beta(:,2) == 1);
Colorednodes{2} = find(beta(:,1) ==1 & beta(:,2) == -1);
Colorednodes{3} = find(beta(:,1) ==-1 & beta(:,2) == 1);
Colorednodes{4} = find(beta(:,1) ==-1 & beta(:,2) == -1);
beta_dist = [1 1; 1 -1; -1 1; -1 -1];

Colorednodes1 = Colorednodes;
Colorednodes{1} = Colorednodes1{2};
Colorednodes{2} = Colorednodes1{1};
beta(Colorednodes{1},:) = 1;
beta(Colorednodes{2},1) = 1;
beta(Colorednodes{2},2) = -1;
beta(Colorednodes{3},2) = 1;
beta(Colorednodes{3},1) = -1;
beta(Colorednodes{4},1) = -1;
beta(Colorednodes{4},2) = -1;