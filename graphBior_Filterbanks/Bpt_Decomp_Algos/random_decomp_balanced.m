function [beta bptG beta_dist Colorednodes] = random_decomp_balanced(A,F)
% A1 = A;
N = length(A);
Fmin = min(F);
F = F - Fmin +1; % minimum color is 1
Fmax = max(F); % maximum value of the color
S = sprintf('%s%d%s','The graph has proper ', Fmax, ' coloring');
disp(S);
nColors = zeros(Fmax,1);
for i = 1:Fmax
    nColors(i) = length(find(F == i));
end
[dontcare I] = sort(nColors,'descend');
tempF = floor((Fmax -4)/4);
if tempF >=0
    for i = 0: tempF
        temp = I(4*i +3);
        I(4*i +3) = I(4*(i+1));
        I(4*(i+1)) = temp;
    end
end
c1 = I(1:2:end);
c2 = I(2:2:end);
Colorednodes = cell(2,1);
Colorednodes{1} = [];
for i = 1:length(c1)
    Colorednodes{1} = union(Colorednodes{1}, find(F == c1(i))); % these nodes have ith color
end
Colorednodes{2} = [];
for i = 1:length(c2)
    Colorednodes{2} = union(Colorednodes{2}, find(F == c2(i))); % these nodes have ith color
end
Fmax =2;
theta = ceil(log2(Fmax)); % this many bpts will be produced
disp('Computing bipartite subgraph decomposition...');
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
