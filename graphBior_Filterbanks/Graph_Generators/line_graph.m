
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
%   Description: This code generates a regular line graph of size N
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A xy] = line_graph(varargin)
switch size(varargin,2)
    case 0 % default config
        N = 50;
    case 1 % size
        N = varargin{1};
    otherwise
        disp('Too many input arguments. generating default case');
        N = 50;
end
% generates a line graph with N nodes
A = diag(ones(N-1,1),-1) + diag(ones(N-1,1),1);
% addpath ..\Vis'ualize Graphs'\
S= sprintf('%s%d','Generating Line Graph with n = ',N);
disp(S)
N = length(A);
% [x,y]= draw_dot(A); % with ovals
theta = 0:1/N:(1-1/N);
theta = theta*2*pi;
% theta =theta(1:100);
% theta = theta*pi;
xy(:,1) = cos(theta);
xy(:,2) = sin(theta);
% figure,
% gplot(G, [x' y'], '.-'); % without ovals
% index = 1:N;
% T = text(x(index)'+0.0001, y(index)'+0.0001,mat2cell(index',ones(1,N),1));