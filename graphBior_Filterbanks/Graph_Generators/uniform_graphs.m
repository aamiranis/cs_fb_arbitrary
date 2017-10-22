
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   "Copyright (c) 2009 The University of Southern California"
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
%   Description: This code generates
%   This code generates a uniform graph from N datapoints uniformly
%   distributed in a [0, 1] x[0,1] field, that contains binary field x.
%   We construct the weighted graph on the nodes
%   based on a thresholded Gaussian kernel weighting function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A xy f] = uniform_graphs(varargin)
if size(varargin,2) == 0 % default
    N = 500;
    sigma = 4;
elseif  size(varargin,2) == 1 % only size is given
    N = cell2mat(varargin(1));
    sigma = 4;
else
    disp('Too many input arguments. Generating default config')
    N = 500;
    sigma = 4;
end

% underlying Field
% Data = imread('sample.png','png');
Data = imread('sample1.jpg','jpeg');
% if (length(size(Data)) == 3)
%     [Data,map] = rgb2ind(Data,10,'dither');
% else
%     [Data,map] = gray2ind(Data,10);
% end
% imshow(Data,map)
% Data = rgb2gray(Data);
Data = double(Data)/255;
[m n] = size(Data);
m = min([m,n]);
Data = Data(1:m,1:m);
n = m;
close all
% Sample N points uniformly in the field
Sample = floor(m*n*rand(1,N))+1; %uniformly sample datapoints
while(length(Sample)~= length(unique(Sample)))
    Sample = floor(m*n*rand(1,N))+1;
end
Sample = Sample';
loc = [floor(Sample/n)+1 , mod(Sample,m)+1];
% figure, imagesc(Data)
% colormap(gray)
% hold on
f = zeros(N,1);
for i = 1:N
    f(i) = Data(loc(i,1),loc(i,2));
    %     plot(loc(i,1),loc(i,2),'ro', 'MarkerSize',10*f(i),  'MarkerFaceColor','r')
end
f = f/norm(f);
% % Compute similarity based on function value
%
% W1 = f*f';
% dW = diag(W1);
% dW = repmat(dW,1,N);
% dW = dW+dW';
% W1 = dW - 2*W1;
% % Dist(Dist < 0.025) = 0;
% % Dist(Dist > 0.2) = 1;
% W1 = exp(-W1)/(2*sigma^2);
% thresh = prctile(W1(:),85);
% W1(W1 <= thresh) = 0;
%
% % compute similarity based on Eucledian embedding
% W2 = loc*loc';
% dW = diag(W2);
% dW = repmat(dW,1,N);
% dW = dW+dW';
% W2 = dW - 2*W2;
% W2 = exp(-W2)/(2*sigma^2);
% thresh = prctile(W2(:),70);
% W2(W2 <= thresh) = 0;
% compute combined similarity based on Eucledian embedding and function
% value
W = [loc f]*[loc f]';
% W = loc*loc';
dW = diag(W);
dW = repmat(dW,1,N);
dW = dW+dW';
W1 = dW - 2*W;
W = exp(-W1);
W = W - diag(diag(W));
thresh = prctile(W(:),85);
W(W <= thresh) = 0;
A = W - diag(diag(W));
A = 0.5*(A+A');
d = sum(A,2);
index = find(~d);
W1 = W1 - diag(diag(W1));
for i = 1:length(index);
    [doncare nearest_node] = sort(W1(index(i),:));
    nearest_node = nearest_node(2);
    A(index(i),nearest_node) =exp(-W1(index(i),nearest_node));
    A(nearest_node,index(i)) = A(index(i),nearest_node);
end
A = double(A>0);
% design adjacency matrix
% A = W2.*W1;



% gplot(A,loc);
addpath ..\Vis'ualize Graphs'\
S= sprintf('%s%d','Generating Unifromly distributed graph with n = ',N);
disp(S)
xy = loc;