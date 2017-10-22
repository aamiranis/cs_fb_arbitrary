
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
%   Description: This code generates a random bipartite graph with
%   paritions of size n and m and edge probability p
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A xy] = random_bipartite_graph(varargin)
switch size(varargin,2)
    case 0 % default config
        m = 50;
        n = 50;
        p = 2*log(n)/(n);
    case 1 % balanced partition graph
        m = varargin{1};
        n = m;
        p = 2*log(n)/(n);
    case 2 %  sizes of two partitions
        m = varargin{1};
        n = varargin{2};
        Max = max([m,n]);
        p = 2*log(Max)/(Max);
    case 3
        m = varargin{1};
        n = varargin{2};
        p = varargin{3};
    otherwise
        disp('Too many input arguments. generating default case');
        m = 50;
        n = 50;
        p = 2*log(n)/(n);
end
G = rand(m,n);
G = double(G < p);
A = zeros(m+n);
A(1:m,m+1:m+n) = G;
A = A + A';
% addpath ..\Vis'ualize Graphs'\
S= sprintf('%s%d%s%d%s%0.5g','Generating Random Bipartite graph with with partition size [n m]= [',n,',',m,'] and edge creation probability p = ',p);
disp(S)
% [x,y]= draw_dot(A); % with ovals
% close
xy = rand(m+n,2);
% xy(:,1) = x(:);
% xy(:,2) = y(:);