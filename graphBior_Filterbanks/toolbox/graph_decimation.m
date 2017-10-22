function [Ad S]= graph_decimation(A,xy,max_level)
% given an adjacency matrix finds maximal independent set based
% decomposition
close all
clc
A = double(A>0); % for now only for unweighted sets
Ad = cell(max_level+1,1);
% xy_d  = cell(max_level+1,1);
S = cell(max_level+1,1);
Ad{1} = A;
S{1} = 1:length(Ad{1});
xy_d{1} = xy;
level = 0;
figure,
gplot(Ad{1},xy_d{1})
V = axis;
S1 = sprintf('%s%d%s', 'Graph at ',level, 'th resolution');
title(S1);
for level = 1:max_level
    N = length(Ad{level});
    % find maximal independent set
    S{level+1} = greedy_indpendent_set(Ad{level});
    Sres = setdiff(1:N,S{level+1}); % remaining nodes
%     tempA =Ad{level}(Sres,S{level+1});
%     d = sum(tempA,2);
%     leftouts = Sres(find(~d));
%     S{level+1} = union(S{level+1},leftouts);
%     Sres = setdiff(Sres,leftouts);
    S{level+1} = Sres;
%     Sres = S{level+1};
    xy_d{level+1} = xy_d{level}(Sres,:);
    % form a bipartite graph
    tempAd = Ad{level};
    tempAd(Sres,Sres) = 0;
    tempAd = tempAd^2;
    tempAd = tempAd - diag(diag(tempAd));
    tempAd = tempAd(Sres,Sres) + Ad{level}(Sres,Sres);
    Ad{level+1} = double(tempAd >0);
    figure,
    gplot(Ad{level+1},xy_d{level+1},'*-')
    S1 = sprintf('%s%d%s', 'Graph at ',level, 'th resolution');
    title(S1);
    axis(V)
end
end

function S = greedy_indpendent_set(A)
% finds maximal independent set
d = sum(A,2);
pendent_nodes = find(~(d-1));
d_max = max(d);
S = [];
while d_max > 0
    d= sum(A,2);
    [d_max loc] = max(d);
    if d_max == 0
        break;
    end
    S = union(S,loc);
    nbrs = find(A(loc,:));
    close_nbrs = union(loc,nbrs);
    A(close_nbrs,:) = 0;
    A(:,close_nbrs) = 0;
end
S = setdiff(S,pendent_nodes); % leave pendent nodes out
end





