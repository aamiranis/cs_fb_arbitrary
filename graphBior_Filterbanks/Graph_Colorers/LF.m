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
%  This file implements a greedy graph coloring algorithm by using 'largest 
%  degree heuristics. Input A is the adjacency matrix. Output F is 
%   the coloring vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = LF(A)
% The Largest Degree First Greedy Graph Coloring Algorithm. Input A is the 
% adjacency matrix. Output F is the coloring vector
N = length(A);
F = zeros(N,1);
deg = sum(A,2);
[dontcare,order] = sort(deg,'descend');

Colors = 0;
for i = 1:length(order)
    node = order(i);
    nbrs = find(A(node,:));
    Rem_colors = setdiff(Colors,union(F(node),F(nbrs)));
    if ~isempty(Rem_colors) % if there is an unused color in existing colors then use it
        F(node) = Rem_colors(1);
    else % if all colors have been used then assign a new color
        max_Color = max(Colors);
        F(node) = max_Color + 1;
        Colors = union(Colors, max_Color +1);
    end
end

        
    
    
