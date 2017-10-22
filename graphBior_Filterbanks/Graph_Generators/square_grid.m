
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
%   Description: This code generates a regular square grid graph with n vertics in each row
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A xy] = square_grid(varargin)
switch size(varargin,2)
    case 0 % default config
        n = 20;
    case 1 % size
        n = varargin{1};
    otherwise
        disp('Too many input arguments. generating default case');
        n = 20;
end
% n = 50;
N = n^2; %total nodes
A = zeros(N);
% interior
for i = 2:n-1
    for j = 2:n-1
        loc1 = (j-1)*n + i;
        loc2 = (j-1)*n + i-1;
        loc3 = (j-1)*n + i+1;
        loc4= (j-2)*n + i;
        loc5 = j*n+i;
        A(loc1,[loc2,loc3,loc4,loc5]) = 1;
    end
    %first column
    loc1 = i;
    loc2 = i-1;
    loc3 = i+1;
    loc4 = n + i;
    A(loc1,[loc2,loc3,loc4]) = 1;
    %last column
    loc1 = (n-1)*n + i;
    loc2 = (n-1)*n + i-1;
    loc3 = (n-1)*n + i+1;
    loc4 = (n-2)*n + i;
    A(loc1,[loc2,loc3,loc4]) = 1;
end
for j = 2:n-1
    % first row
    loc1 = (j-1)*n + 1;
    loc2 = (j-2)*n + 1;
    loc3 = j*n + 1;
    loc4 = (j-1)*n + 2;
    A(loc1,[loc2,loc3,loc4]) = 1;
    %last row
    loc1 = (j-1)*n + n;
    loc2 = (j-2)*n + n;
    loc3 = j*n + n;
    loc4 = (j-1)*n + n-1;
    A(loc1,[loc2,loc3,loc4]) = 1;
end

% corner points
% top left
loc1 = 1;
loc2 = 2;
loc3 = n+1;
A(loc1,[loc2,loc3]) = 1;
% top right
loc1 = (n-1)*n + 1;
loc2 = (n-2)*n + 1;
loc3 = (n-1)*n +2;
A(loc1,[loc2,loc3]) = 1;
%bottom left
loc1 =  n;
loc2 = n-1;
loc3 = 2*n;
A(loc1,[loc2,loc3]) = 1;
%bottom right
loc1 =  n^2;
loc2 = n^2 -1;
loc3 = n*(n-1);
A(loc1,[loc2,loc3]) = 1;
% G = G + G';
% G = double(G>0);
addpath ..\Vis'ualize Graphs'\
S= sprintf('%s%d','Generating Square Grid with side length n = ',n);
disp(S)
% [x,y]= draw_dot(A); % with ovals
x = 1:n;
x = repmat(x(:),1,n);
y = 1:n;
y = repmat(y(:)', n ,1);
xy(:,1) = x(:);
xy(:,2) = y(:);
% figure,
% gplot(G, [x' y'], '.-'); % without ovals
% index = 1:N;
% T = text(x(index)'+0.0001, y(index)'+0.0001,mat2cell(index',ones(1,N),1));
