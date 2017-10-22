
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
%   Description: This code generates an Erdos Renyi uniform random
%   graph with n vertics and edge probability p
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A xy] = Erdos_random_graph(varargin)
if size(varargin,2) == 0 % default
    n = 50;
    p = 0.15;
elseif  size(varargin,2) == 1 % only size is given
    n = varargin{1};
    p = 2*log(n)/n;
elseif size(varargin,2) == 2
    n = varargin{1};
    p = varargin{2};
else
    disp('Too many input arguments. Generating default config')
    n = 50;
    p = 0.15;    
end


G = rand(n);
G = G < p;
G = triu(G,1);
G = G + G';
Sum = sum(G,2);
ind = find(Sum ~=0);
G = G(ind,ind);
% G = connected_component(G);
n = length(G);
addpath ..\Vis'ualize Graphs'\
S= sprintf('%s%d%s%0.5g','Generating Erdos-Reyni Random Graph with n = ',n,' and p = ',p);
disp(S)
[x,y]= draw_dot(G); % with ovals
A = G;
xy(:,1) = x(:);
xy(:,2) = y(:);