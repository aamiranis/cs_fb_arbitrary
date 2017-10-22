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
%  This file implements a greedy graph coloring algorithm by using 'degree 
%  of saturation heuristics. Input A is the adjacency matrix. Output F is 
%   the coloring vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = DSATUR(A)
N = length(A);
F = zeros(N,1);
deg = sum(A,2);
[dontcare,node] = max(deg); 
F(node) = 1;
Colors = [0 1];
for i = 2:N % color remaining nodes
    uncolored_nodes = find(~F);
    colored_nodes = find(F);
    deg_s = sum(A(uncolored_nodes,colored_nodes),2);
    [dontcare,loc] = max(deg_s);
    node = uncolored_nodes(loc(1));
%     if node == 2418 || node == 2389
%         ghdgh = 0;
%     end
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

        
    
    
