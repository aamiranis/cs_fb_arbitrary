
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
%   Description: This code simulates a scale free small world model based on 
%  " Barabasi-Albert (BA) model". The data generator constructs a graph 
%   G = (V,E) as follows: The graph begins with an initial network of 
%   m0>=2 nodes. It should be noted that  and the degree of each node in the 
%   initial network should be at least 1, otherwise it will always 
%   remain disconnected from the rest of the network.
%   New nodes are added to the network one at a time. Each new node 
%   is connected to m of the existing with a probability that is biased 
%   so that it is proportional to the number of links that the existing
%   node already has. Formally, the probability pi that the new node is connected to node 
%   i is ki/sum(kj) ,where ki is the degree of node i. Heavily linked nodes (
%   "hubs") tend to quickly accumulate even more links, while nodes with only 
%   a few links are unlikely to be chosen as the destination for a new link. 
%   The new nodes have a "preference" to attach themselves to the already 
%  heavily linked nodes.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A xy] = small_world_BA(varargin)
switch size(varargin,2)
    case 0 % default config
        N = 100;
        m0 = 10;
    case 1 % size 
        N = varargin{1};
        m0 = 10;
    case 2 %  sizes 
        N = varargin{1}; % total number of nodes
        m0 = varargin{2}; % initial size of the network
    otherwise
        disp('Too many input arguments. generating default case');
        N = 100;
        m0 = 10;
end

% INPUTS
%N = # of nodes
epsilon = 10^-3;
adj = zeros(N);
for i = 1:m0-1
    adj(i,i+1) = 1;
    adj(i+1,i) = 1;
end
for i = m0+1:N
    sumk = sum(sum(adj(1:i-1,1:i-1)))/2;
    P = sum(adj(1:i-1,:),2)/sumk;
    r = rand(i-1,1);
    index = r <= P;
    adj(i,index) = 1;
    adj(index,i) = 1;
end 

% now take out the largest connected component out of it
index = find(sum(adj,2)>0);
adj = adj(index,index); 

addpath ..\Vis'ualize Graphs'\
S= sprintf('%s%d','Generating random BA graph with number of nodes = ', N);
disp(S)
[x,y,names]= draw_dot(adj); % with ovals
A = adj;
xy(:,1) = x(:);
xy(:,2) = y(:);
% save ..\T'ree wavelets on graphs ver 1'\smallworld_data  A x y labels class_sizes


