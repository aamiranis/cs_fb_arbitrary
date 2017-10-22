function [beta bptG beta_dist Colorednodes] = harary_decomp1(A,F)
% A1 = A;
N = length(A);
Fmin = min(F);
F = F - Fmin +1; % minimum color is 1
Fmax = max(F); % maximum value of the color
S = sprintf('%s%d%s','The graph has proper ', Fmax, ' coloring');
disp(S);
theta = ceil(log2(Fmax)); % this many bpts will be produced
disp([num2str(theta),' bipartite subgraphs will be generated.']);
disp('Computing bipartite subgraph decomposition...');
new_order =  dec2bin(0:(2^(theta)-1));
new_order = fliplr(new_order);
new_order = bin2dec(new_order);
new_order = new_order(new_order < Fmax) + 1;
Colorednodes = cell(Fmax,1);
for i = 1:Fmax
    Colorednodes{i} = find(F == new_order(i)); % these nodes have ith color
end

beta_dist = dec2bin(sort_colors(Fmax,theta));
beta_dist = double(beta_dist) - double('0');
beta_dist = 2*beta_dist - 1;
beta = zeros(N,theta);
bptG = zeros(N,N,theta);
for i = 1:theta
    for j = 1: Fmax
        if beta_dist(j,i) == 1
            beta(Colorednodes{j},i) = 1;
        else
            beta(Colorednodes{j},i) = -1;
        end
    end
    S1 = find(beta(:,i) == 1);
    S2 = find(beta(:,i) == -1);
    bptG(S1,S2,i) = A(S1,S2);
    bptG(S2,S1,i) = A(S2,S1);
    A(S1,S2) = 0;
    A(S2,S1) = 0;
end
edges = sum(sum(bptG,1),2)/2;
edges = edges(:);
cum_edges = cumsum(edges)/sum(edges);
theta1 = find(cum_edges > 0.90, 1,'first');
theta = theta1;
% theta =2;
bptG = bptG(:,:,1:theta);
beta = beta(:,1:theta);
beta_dist = beta_dist(:,1:theta);
beta_dist1 = 0.5*(beta_dist + 1);
b = 2.^((theta-1):-1:0);
b = b';
Fmax1 = 2^theta;
F = Fmax1 - beta_dist1*b;
Colorednodes1 = cell(Fmax1,1);
for i = 1: Fmax1
    Colorednodes1{i} = [];
end
for i = 1: Fmax
    Colorednodes1{F(i)} = union(Colorednodes1{F(i)}, Colorednodes{i});
end
Colorednodes = Colorednodes1;
beta_dist = dec2bin(sort_colors(Fmax1,theta));
beta_dist = double(beta_dist) - double('0');
beta_dist = 2*beta_dist - 1;
    

