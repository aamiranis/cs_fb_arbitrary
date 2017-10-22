
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
function [bptG Colorednodes beta_dist loc] = image_graphs_multi_Bior(Data,max_level,edgemap,norm_type,filterlen)
link_weight_low = 1/10;
theta = 2;
s_im = length(Data);
N = s_im*s_im;
bptG = cell(max_level,2);
Colorednodes = cell(max_level,4);
beta_dist = cell(max_level,1);
loc = cell(max_level,1);
beta(:,1) = image_downsampling_fn(s_im,'rectangle');
beta(:,2) = image_downsampling_fn(s_im,'diamond');
beta_temp = 0.5*(1+beta);
F = 2*beta_temp(:,1 ) + beta_temp(:,2);
clear beta_temp
F = 2^theta - F ;
Fmax = 2^theta;
for i = 1:Fmax
    Colorednodes{1,i} = find(F == i); % these nodes have ith color
end
beta_dist{1} = double(dec2bin(0:Fmax-1)) - double('0');
beta_dist{1} = 1 - 2*beta_dist{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First create first level image-graph without using the edge-map
dim = [s_im s_im];
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

% % A1 = exp(-1)*A1;
% tempD1 = (f_bw(cind1) - f_bw(nind1)).^2;
% tempD1 = exp(-tempD1);
% % A1 = exp(-1)*A1;
% tempD2 = (f_bw(cind2) - f_bw(nind2)).^2;
% tempD2 = exp(-tempD2);
% 
% thresh = 0.5;%prctile(tempD(:),40);
% tempD1(tempD1 <= thresh) = 0;
% tempD2(tempD2 <= thresh) = 0;

A1=sparse([cind1,nind1],[nind1,cind1],[ones(numel(cind1),1);ones(numel(cind1),1)],N,N);
A2=sparse([cind2,nind2],[nind2,cind2],[ones(numel(cind2),1);ones(numel(cind2),1)],N,N);



%% Implemenent undecimated HH band filter 
%% Section 2: Filterbank implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Normalized Laplacian Matrices for Each Bpt graph
% disp('Computing normalized Laplacian matrices for each subgraph...');
switch norm_type
    case 'sym'
        d1 = sum(A1,2);
        d1(d1 == 0) = 1; % for isolated nodes
        d1_inv = d1.^(-0.5);
        D1_inv = spdiags(d1_inv, 0, N, N);
        An = D1_inv*A1*D1_inv;
        An = 0.5*(An + An');
        Ln1 = speye(N) - An;
        d1 = sum(A2,2);
        d1(d1 == 0) = 1; % for isolated nodes
        d1_inv = d1.^(-0.5);
        D1_inv = spdiags(d1_inv, 0, N, N);
        An = D1_inv*A2*D1_inv;
        An = 0.5*(An + An');
        Ln2 = speye(N) - An;
    case 'asym' 
        d1 = sum(A1,2);
        d1(d1 == 0) = 1; % for isolated nodes
        d1_inv = d1.^(-1);
        D1_inv = spdiags(d1_inv, 0, N, N);
        An = D1_inv*(0.5*(A1 + A1'));
        Ln1 = speye(N) - An;
        d1 = sum(A2,2);
        d1(d1 == 0) = 1; % for isolated nodes
        d1_inv = d1.^(-1);
        D1_inv = spdiags(d1_inv, 0, N, N);
        An = D1_inv*(0.5*(A2 + A2'));
        Ln2 = speye(N) - An;
    otherwise
        disp('Unknown normalization option')
        return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% design biorthogonal kernels
N_lo = floor(filterlen/2) + mod(filterlen,2);
N_hi = N_lo;
% S = sprintf('%s%d%s ' , ' Computing a ',filterlen,' ^th order approximation of Meyer kernel');
% disp(S)
[lo_d,hi_d] = biorth_kernel(N_lo,N_hi); % Nd zeros for lowpass Nr zeros for highpass
filterlen_hi = length(roots(hi_d));
filterlen_lo = length(roots(lo_d));
h0 = @(x)(polyval(lo_d,x));
h1 = @(x)(polyval(hi_d,x));
g0 = @(x)(polyval(hi_d,2 - x));
g1 = @(x)(polyval(lo_d,2 - x));
arange = [0 2];
c_d{1}=sgwt_cheby_coeff(h0,filterlen_lo,filterlen_lo+1,arange);
c_d{2}=sgwt_cheby_coeff(h1,filterlen_hi,filterlen_hi+1,arange);
c_r{1}=sgwt_cheby_coeff(g0,filterlen_hi,filterlen_hi+1,arange);
c_r{2}=sgwt_cheby_coeff(g1,filterlen_lo,filterlen_lo+1,arange);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Filterbank Output at each channel
% Compute Filterbank Output at each channel
f = double(Data(:));
f = f/norm(f);
f_w = zeros(N,4);
tempf_w = f;
tempf_w= sgwt_cheby_op(tempf_w,Ln1,c_d{1},arange);
tempf_w = sgwt_cheby_op(tempf_w,Ln2,c_d{1},arange);
f_w(Colorednodes{1,1},1) = tempf_w(Colorednodes{1,1});
tempf_w = f;
tempf_w= sgwt_cheby_op(tempf_w,Ln1,c_d{1},arange);
tempf_w = sgwt_cheby_op(tempf_w,Ln2,c_d{2},arange);
f_w(Colorednodes{1,2},1) = tempf_w(Colorednodes{1,2});
tempf_w = f;
tempf_w= sgwt_cheby_op(tempf_w,Ln1,c_d{2},arange);
tempf_w = sgwt_cheby_op(tempf_w,Ln2,c_d{1},arange);
f_w(Colorednodes{1,3},1) = tempf_w(Colorednodes{1,3});
tempf_w = f;
tempf_w= sgwt_cheby_op(tempf_w,Ln1,c_d{2},arange);
tempf_w = sgwt_cheby_op(tempf_w,Ln2,c_d{2},arange);
f_w(Colorednodes{1,4},1) = tempf_w(Colorednodes{1,4});

% Reconstruct low-pass and high-pass 
tempf_hat = f_w(:,1);
tempf_hat = sgwt_cheby_op(tempf_hat,Ln2,c_r{1},arange);
f_low = sgwt_cheby_op(tempf_hat,Ln1,c_r{1},arange);
f_high = f - f_low;
BW = f_high;
BW = reshape(abs(BW),s_im,s_im);
[histN, X] = hist(BW(:),256);
histN = cumsum(histN)/sum(histN);
pos = find(histN > (0.6 + 0.4*histN(1)),1);
BW = im2bw(BW,X(pos));
% % Create Edge map from the image
% BW = edge(Data,'log');
cc = bwconncomp(BW);
numPixels = cellfun(@numel,cc.PixelIdxList);
small_cc = find(numPixels < 50);
for i = 1:length(small_cc)
    idx = small_cc(i);
    BW(cc.PixelIdxList{idx}) = 0;
end
BW = double(BW);
SE = strel('square',2);
BW = imdilate(BW,SE);
BW = imerode(BW,SE);
BW = imerode(BW,SE);
BW = imerode(BW,SE);
BW = imerode(BW,SE);
BW1 = BW;
if edgemap == 0
    BW = 0*BW;
end
save BWdata BW1
N = s_im*s_im;
index = repmat(1:s_im,[s_im 1]);
loc{1} = zeros(N,2);
loc{1}(:,1) = index(:);
index = index';
loc{1}(:,2) = index(:)';
close all
f_bw = double(BW(:));
f = Data(:);
% figure,
% imagesc(BW);
% colormap(gray);
% f = f/norm(f);


% rectangular connectivity of edge pixels
thresh = 0.5;
[edge_i edge_j] = find(BW);
ci = edge_i;
cj = edge_j;
ni = [edge_i (edge_i + 1)  (edge_i -1) edge_i];
nj = [(edge_j+1) edge_j  edge_j  (edge_j -1)];
cind3= [];
nind3 =[];
for i = 1:length(ci)
    valid=(ni(i,:)>=1 & ni(i,:)<=dim(1) & nj(i,:)>=1 & nj(i,:)<=dim(2));
    temp_cind = dim(1)*(cj(i)-1)+ci(i);
    temp_nind = dim(1)*(nj(i,valid)-1)+ni(i,valid);
    %     nind3=[nind3 ; temp_nind(:)];
    
    %     cind3=[cind3 ; temp_cind(:)];
    temp_w = (f_low(temp_cind) - f_low(temp_nind)).^2;
    temp_w = temp_w/(sum(temp_w) + 10^-7);
    [temp_w order] = sort(temp_w);
    [diffw_mag diffw] = max(abs(diff(temp_w)));
    if diffw_mag > 10^-2
        temp_nind = temp_nind(order((diffw+1):end));
        temp_cind = repmat(temp_cind,numel(temp_nind),1);
        cind3 = [cind3; temp_cind(:)];
        nind3 = [nind3; temp_nind(:)];
    end
    %     tempD3 = [tempD3 ; temp_w(:)];
end
% A1(cind,nind) = link_weight_low;
% A1(nind,cind) = link_weight_low;
A3=sparse([cind3; nind3],[nind3;cind3],ones(2*numel(cind3),1),N,N);
A3 = (-1+link_weight_low)*double(A3 >0);
A1 = A1+A3;
% A3=sparse([cind3,nind3],[nind3,cind3],[tempD3(:);tempD3(:)],N,N);
% index = [cind3 nind3; nind3 cind3];
% tempD3 = [tempD3(:);tempD3(:)];
% [index I J] = unique(index,'rows');
% tempD3 = tempD3(I);
% A3=sparse(index(:,1),index(:,2),tempD3(:),N,N);
% A3 = triu(A3);
% A3 = A3+A3';

% diamond connectivity of edge pixels
[edge_i edge_j] = find(BW);
ci = edge_i;
cj = edge_j;
ni = [(edge_i+1) (edge_i -1) (edge_i +1) (edge_i -1)];
nj = [(edge_j +1) (edge_j -1) (edge_j -1) (edge_j+1)];
cind4= [];
nind4 =[];
for i = 1:length(ci)
    valid=(ni(i,:)>=1 & ni(i,:)<=dim(1) & nj(i,:)>=1 & nj(i,:)<=dim(2));
    temp_nind = dim(1)*(nj(i,valid)-1)+ni(i,valid);
    temp_cind = dim(1)*(cj(i)-1)+ci(i);
    temp_w = (f(temp_cind) - f(temp_nind)).^2;
    temp_w = temp_w/(sum(temp_w) + 10^-7);
    [temp_w order] = sort(temp_w);
    [diffw_mag diffw] = max(abs(diff(temp_w)));
    if diffw_mag > 10^-2
        temp_nind = temp_nind(order((diffw+1):end));
        temp_cind = repmat(temp_cind,numel(temp_nind),1);
        cind4 = [cind4; temp_cind(:)];
        nind4 = [nind4; temp_nind(:)];
    end
end
A4=sparse([cind4; nind4],[nind4;cind4],ones(2*numel(cind4),1),N,N);
A4 = (-1+link_weight_low)*double(A4 >0);
A2 = A2+A4;
% A1 = A1 - A1 .*( A3 >0) + A3;
% A1 = double(A1>0);
bptG{1,1} = 0.5*(A1+A1');
% A2 = A2 - A2.* ( A4>0) + A4;
% A2 = double(A2 >0);
bptG{1,2} = 0.5*(A2+A2');


%% Create image-graphs in the subsequent levels using the original graphs.
for level = 2:max_level
    LL = Colorednodes{level-1,1};
    N = length(LL);
    s_im = sqrt(N);
    bptG{level,1} = bptG{level-1,1}^2;
    bptG{level,1} = bptG{level,1}(LL,LL);
    bptG{level,1} = bptG{level,1} - spdiags(spdiags(bptG{level,1},0),0,N,N);
    %     D = sum(bptG{level,1},2);
    bptG{level,2} = bptG{level-1,2}^2;
    bptG{level,2} = bptG{level,2}(LL,LL);
    bptG{level,2} = bptG{level,2} - spdiags(spdiags(bptG{level,2},0),0,N,N);
%     bptG{level,2} = double(bptG{level,2} >0);
    %%% remove rectangular links from bptG
    % rectanglar links
    dim = [s_im s_im];
    [alli,allj]=find(ones(dim));
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
    tempA=sparse([cind1,nind1],[nind1,cind1],ones(1,2*numel(ni)),N,N);
    tempA = bptG{level,2}.*tempA; % rectangular connectivity
%     bptG{level,1} = bptG{level,1} + tempA;
%     bptG{level,1} = double(bptG{level,1} >0);
    bptG{level,2} = bptG{level,2} -tempA; % keep only diamond connections
    loc{level} = loc{level-1}(LL,:);
    clear beta
    beta(:,1) = image_downsampling_fn(s_im,'rectangle');
    beta(:,2) = image_downsampling_fn(s_im,'diamond');
    beta_temp = 0.5*(1+beta);
    F = 2*beta_temp(:,1 ) + beta_temp(:,2);
    clear beta_temp
    F = 2^theta - F ;
    Fmax = 2^theta;
    for i = 1:Fmax
        Colorednodes{level,i} = find(F == i); % these nodes have ith color
    end
    beta_dist{level} = double(dec2bin(0:Fmax-1)) - double('0');
    beta_dist{level} = 1 - 2*beta_dist{level};
     
%     %% treat isolated pixels in the diamond graph here
%     S3 = union(Colorednodes{level,1},Colorednodes{level,3});
%     S4 = union(Colorednodes{level,2},Colorednodes{level,4});    
%     tempd =sum(bptG{level,2},2);
%     cind = find(~tempd);
%     cind1 = intersect(cind,S3);
%     cind2 = intersect(cind,S4);
%     [row col] =  find(tempA);
%     for i = 1:length(cind1)
%         temploc = ~(row -cind1(i));
%         nind = col(temploc);
%         %         choice1 = nind(1);
%         nind = setdiff(nind,S3);
%         if ~isempty(nind)
%             bptG{level,2}(cind1(i),nind(1)) =  1;
%             bptG{level,2}(nind(1),cind1(i)) =  1;
%         end
%     end
%     for i = 1:length(cind2)
%         temploc = ~(row -cind2(i));
%         nind = col(temploc);
%         %         choice1 = nind(1);
%         nind = setdiff(nind,S4);
%         if ~isempty(nind)
%             bptG{level,2}(cind2(i),nind(1)) =  1;
%             bptG{level,2}(nind(1),cind2(i)) =  1;
%         end
%     end
%     tempd =sum(bptG{level,1},2);
%     cind = find(~tempd);
%     cind = setdiff(cind, intersect(cind,Colorednodes{level,1}));
%     isolated{level} =cind;
%     tempd =sum(bptG{level,2},2);
%     cind = find(~tempd);
%     cind = setdiff(cind, intersect(cind,Colorednodes{level,1}));
%     isolated{level} =union(isolated{level},cind);
end
% scrsz = get(0,'ScreenSize');
% height = scrsz(4)/2;
% width =  scrsz(3)/2;
% xinit = 30;
% yinit = 30;
% for level = 1:max_level
%     figure, gplot(bptG{level,1},loc{level});
%     hold on
%     gplot(A3,loc{level},'r-');
%     set(gca,'YDir','reverse');
%     figure, gplot(bptG{level,2},loc{level});
%     hold on
%     gplot(A4,loc{level},'r-');
%     set(gcf,'Position',[xinit,yinit,width,height]);
%     set(gca,'YDir','reverse')
%     height = height/2;
%     width =  width/2;
%     xinit = xinit +width;
%     yinit = yinit+height;
% end



%%

