
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
function [A1 A2 xy f BW] = image_graphs(filename,filetype)
% underlying Field
% Data = imread('sample1.jpg','jpeg');
% Data = imread('sample1.png','png');
% filename = 'lena.jpg';
Data = imread(filename,filetype);
if length(size(Data)) == 3
    Data = rgb2gray(Data);
end
Data = double(Data);%/255;
[m n] = size(Data);
m = min([m,n]);
m = floor(m/2)*2;
Data = Data(1:m,1:m);
BW = edge(Data,'prewitt');
BW = double(BW);
n = m;
N = m*m;
index = repmat(1:m,[m 1]);
loc = zeros(N,2);
loc(:,1) = index(:);
index = index';
loc(:,2) = index(:)';
close all
f_bw = double(BW(:));
f = Data(:);
% figure,
% imagesc(BW);
% colormap(gray);
% f = f/norm(f);
xy = loc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% copied from sgwt meshmat
dim = [m m];
[alli,allj]=find(ones(dim));

% rectanglar links
ci=[alli;alli];
cj=[allj;allj];
ni=[alli  ; alli+1];
nj=[allj+1; allj];
% prune edges at boundary
valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2));
ni=ni(valid);
nj=nj(valid);
ci=ci(valid);
cj=cj(valid);
cind1=dim(1)*(cj-1)+ci;
nind1=dim(1)*(nj-1)+ni;

% diagonal links
ci=[alli;alli];
cj=[allj;allj];
ni=[alli+1  ; alli+1];
nj=[allj-1; allj+1];
% prune edges at boundary
valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2));
ni=ni(valid);
nj=nj(valid);
ci=ci(valid);
cj=cj(valid);
cind2=dim(1)*(cj-1)+ci;
nind2=dim(1)*(nj-1)+ni;

% A1 = exp(-1)*A1;
tempD1 = (f_bw(cind1) - f_bw(nind1)).^2;
tempD1 = exp(-tempD1);
% A1 = exp(-1)*A1;
tempD2 = (f_bw(cind2) - f_bw(nind2)).^2;
tempD2 = exp(-tempD2);

thresh = 0.5;%prctile(tempD(:),40);
tempD1(tempD1 <= thresh) = 0;
tempD2(tempD2 <= thresh) = 0;

A1=sparse([cind1,nind1],[nind1,cind1],[tempD1(:);tempD1(:)],N,N);
A2=sparse([cind2,nind2],[nind2,cind2],[tempD2(:);tempD2(:)],N,N);


%% Edge Pixels are treated here
% rectangular connectivity of edge pixels
[edge_i edge_j] = find(BW);
ci = edge_i;
cj = edge_j;
ni = [edge_i (edge_i + 1)  (edge_i -1) edge_i];
nj = [(edge_j+1) edge_j  edge_j  (edge_j -1)];
cind3= [];
nind3 =[];
tempD3 = [];
for i = 1:length(ci)
    valid=(ni(i,:)>=1 & ni(i,:)<=dim(1) & nj(i,:)>=1 & nj(i,:)<=dim(2));
    temp_nind = dim(1)*(nj(i,valid)-1)+ni(i,valid);
    nind3=[nind3 ; temp_nind(:)];
    v = length(temp_nind);
    temp_cind = repmat(dim(1)*(cj(i)-1)+ci(i),v,1);
    cind3=[cind3 ; temp_cind(:)];
    temp_w = (f(temp_cind) - f(temp_nind)).^2;
    if max(temp_w)>10^-10
        temp_w = temp_w/max(temp_w);
    end
    tempD3 = [tempD3 ; temp_w(:)];
end
thresh = 0.9;
tempD3(tempD3 <= thresh) = 0;
A3=sparse([cind3,nind3],[nind3,cind3],[tempD3(:);tempD3(:)],N,N);
    
% diamond connectivity of edge pixels
[edge_i edge_j] = find(BW);
ci = edge_i;
cj = edge_j;
ni = [(edge_i+1) (edge_i -1) (edge_i +1) (edge_i -1)];
nj = [(edge_j +1) (edge_j -1) (edge_j -1) (edge_j+1)];
cind4= [];
nind4 =[];
tempD4 = [];
for i = 1:length(ci)
    valid=(ni(i,:)>=1 & ni(i,:)<=dim(1) & nj(i,:)>=1 & nj(i,:)<=dim(2));
    temp_nind = dim(1)*(nj(i,valid)-1)+ni(i,valid);
    nind4=[nind4 ; temp_nind(:)];
    v = length(temp_nind);
    temp_cind = repmat(dim(1)*(cj(i)-1)+ci(i),v,1);
    cind4=[cind4 ; temp_cind(:)];
    temp_w = (f(temp_cind) - f(temp_nind)).^2;
    if max(temp_w)>10^-10
        temp_w = temp_w/max(temp_w);
    end
    tempD4 = [tempD4 ; temp_w(:)];
end
thresh = 0.9;
tempD4(tempD4 <= thresh) = 0;
A4=sparse([cind4,nind4],[nind4,cind4],[tempD4(:);tempD4(:)],N,N);
    
A1 = A1 + A3;
A1 = double(A1>0);
A1 = 0.5*(A1+A1');
A2 = A2 + A4;
A2 = double(A2 >0);
A2 = 0.5*(A2+A2');

