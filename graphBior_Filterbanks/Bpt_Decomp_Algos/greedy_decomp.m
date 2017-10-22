function [beta bptG beta_dist Colorednodes] = greedy_decomp(A)

[m n] = size(A);
% randomly assign 4 coloring
label = 4*rand(n,1);
label = floor(label);
label(find(label ==4)) = 3;
label = label+1;
% label = repmat(label',n,1);
% L = label.*A;
% conflict = (L == label');
% conflict = conflict.*A;
% C(1) = sum(sum(conflict))/sum(sum(A));
% iterations
max_iter = 40;
for iter = 1: max_iter
    active = rand(n,1);
    active = active <0.3; % activate 30% of the node at a time
    F = find(active == 1);
    for i = 1:length(F)
        k = F(i);
        temp1 = sum(conflict(k,:)); % number of conflicts
        temp2 = sum(A(k,:)); % number of nbr edges
        if temp1/temp2 > 0.5 % if conflicts are more
            label(:,k) = -label(:,k);
            L(:,k) = label(:,k).*A(:,k);
            conflict = (L == label');
            conflict = conflict.*A;
        end
    end
    C(iter) = sum(sum(conflict))/sum(sum(A));
end
label = label(1,:)';
conflict = C(max_iter)*sum(sum(A));
conflict = conflict/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%