%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   "Copyright (c) 2011 The University of Southern California"
%   All rights reserved.
%
%   Permission to use, copy, modify, and distribute this software and its
%   documentation for any purpose, without fee, and without written
%   agreement is hereby granted, provided that the above copyright notice,
%   the following two paragraphs and the author appear in all copies of
%   this software.
%
%   NO REPRESENTATIONS ARE MADE ABOUT THE SUITABILITY OF THE SOFTWARE
%   FOR ANY	PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED
%   WARRANTY.
%
%   Neither the software developers, the Compression Research Group,
%   or USC, shall be liable for any damages suffered from using this
%   software.
%
%   Author: Sunil K Narang
%   Director: Prof. Antonio Ortega
%   Compression Research Group, University of Southern California
%   http://biron.usc.edu/wiki/index.php?title=CompressionGroup
%   Contact: kumarsun@usc.edu
%
%   Date last modified:	07/05/2011 kumarsun
%
%   Description:
% This function implements Backtracing sequential coloring described by 
% Walter Klotz in the paper:  "Graph Coloring Algorithms"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = BSC(A)
k = 1;
N = length(A)
F = zeros(1,N);
deg = sum(A,2);
degs= 0*deg;
[deg, order]= sort(deg,'descend');
start = 1;   % starting index
optColorNumber = deg(1)+1; % since all graphs can be colored with (max degree +1) colors
x = order(1); % current vertex to load
colors(1) = 0;
U = 1;  % variable for the set of free colors
freeColors = cell(1,N);
freeColors{x} = U; % set of free colors of x

while(start >= 1)
    back = 0;
    for i = start:N
        if i> start
            x = find(~F);
            if isempty(x)
                break;
            end
            y = find(degs(x) == max(degs(x)));
            x = x(y(1));
            nbrs = find(A(x,:));
            forbidden_colors = unique(F(nbrs));
            U = setdiff(1:optColorNumber, forbidden_colors);
            freeColors{x} = U;
        end
        if ~isempty(U) > 0
            s = U(1);
            F(x) = s;
            nbrs = find(A(x,:));
            for j = 1: length(nbrs)
                nbrs2 = find(A(nbrs(j),:));
                cols = unique(F(nbrs2));
                cols = setdiff(cols,0);
                degs(nbrs(j)) = length(cols);
            end
            U = setdiff(U,s);
            freeColors{x} = U;
            l = 0;
            if i > 1
                l = colors(i-1);
            end
            colors(i) = max(s,l);
        else
            start = i -1;
            back = 1;
            break;
        end
    end
    if true(back)
        if start >=1
            x = order(start);
            F(x) = 0;
            U = freeColors{x};
        end
    else
        Fopt = F;
        optColorNumber = colors(N-1);
        i = find(F == optColorNumber);
        if ~isempty(i)
            i = i(1);
        else
            i = 0;
        end
        start = i-1;
        if start < 1
            break;
        end
        F(start:N) = 0;
        for i = 1 :start
            x = order(i);
            U = freeColors{x};
            tempU = find(U >= optColorNumber);
            U = setdiff(U,tempU);
            freeColors{x} = U;
        end
    end
end
F = Fopt;










