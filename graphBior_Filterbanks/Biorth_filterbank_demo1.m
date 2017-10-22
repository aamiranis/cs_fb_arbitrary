%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   "Copyright (c) 2012 The University of Southern California"
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
%   Date last modified:	02/01/2012 Sunil K. Narang
%
%   Description:
%   This file generates wavelet coefficients using
%   graphBior filterbanks on the 8-connected graph formulation
%   of 2D digital images as proposed in the paper:
%   "S. K. Narang and Antonio Ortega, "Compact Support Biorthogonal Wavelet Filterbanks
%    for Arbitrary Undirected Graphs" available as preprint version arXiv:1210.8129,
%%
%   Uses:
%   W = Biorth_filterbank_forward(); % compute wavelet coefficients on default image
%   W = Biorth_filterbank_forward('filename'); % compute wavelet coefficients on image <filename>
%   W stores the wavelet coefficients in a pyramid layout.
%
%   Optional parameter-value pairs:
%
%   maxLevel: maximum number of decomposition levels (default 3)
%   filterLen: length of graphBior filters (default 8)
%   Ltype: type of Laplacian matrix used 'sym'/'asym' (default 'asym')
%   GC : enable gain compensation module 'on'/'off' (default 'on')
%   edgemap: edge-aware transform 0/1 (default 1)
%   nnz_factor: fraction of non-zero highpass coefficients used for
%   reconstruction (all lowpass coefficients at the coarsest resolution are
%   used)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [wav_coeffs] = Biorth_filterbank_demo1(varargin)
function [] = Biorth_filterbank_demo1(varargin)
clc
close all
disp('***********Welcome to Demo 1 of Graph-QMF*********')
disp('In this demo we implement a 2-dim separable graph-QMF filterbank on any Image')

addpath Graph_Generators/
addpath Graph_kernels/
addpath sgwt_toolbox/
addpath toolbox/
addpath exportfig/
addpath Datasets/
p = inputParser;
p.addOptional('filename','coins.png',@(f) exist(f,'file'));
p.addParamValue('maxLevel',4,@isscalar);
p.addParamValue('filterLen',10,@isscalar);
p.addParamValue('LType','asym',@isstr);
p.addParamValue('GC','on',@isstr);
p.addParamValue('edgemap',1,@islogical);
p.addParamValue('nnz_factor',0.03,@isreal);
p.parse(varargin{:});
opt = p.Results;

% Set parameters
filename = opt.filename;
tempI = imfinfo(filename);
filetype = tempI.Format;
max_level = opt.maxLevel; % number of decomposition levels
filterlen = opt.filterLen; % filter length
norm_type = opt.LType;
GC_module = opt.GC;
edgemap = opt.edgemap;
nnz_factor = opt.nnz_factor; % fraction of non-zero highpass coefficients
if isempty(nnz_factor)
    nnz_factor = 1;%(4^(max_level-2));
end
if nnz_factor <0 
    nnz_factor = 0;
end
if nnz_factor >1 
    nnz_factor = 1;
end
    
S = sprintf('%s%d\n%s%d\n%s%s\n%s%d\n%s%s', '# decomposition levels = ', max_level,'Filterlength = ', filterlen,'GC module is ',GC_module, 'Using edgemap = ', edgemap,'Type of Laplacian = ',norm_type);
disp(S);

theta = 2; % number of bipartite graphs
Fmax = 2^theta; % Graph Coloring
N = zeros(max_level,1); %size of bipartite graphs at each level

%% Section 1: Image Graph Formulation
% Graph Signal
Data = imread(filename,filetype);
if length(size(Data)) == 3
    Data = rgb2gray(Data);
end
[m n] = size(Data);
m = min([m,n]);
s_im = floor(m/2^max_level)*2^max_level;
% s_im = 256;
Data = Data(1:s_im,1:s_im);
f = double(Data(:));
% +rand(s_im^2,1);
% f = f/norm(f);
% Graphs
[bptG Colorednodes, beta_dist loc] = image_graphs_multi(Data,max_level,edgemap);  % image graphs


%% Section 2: Filterbank implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Normalized Laplacian Matrices for Each Bpt graph
disp('Computing normalized Laplacian matrices for each subgraph...');
Ln_bpt = cell(max_level,theta);
switch norm_type
    case 'sym'
        for level = 1:max_level
            N(level) = length(bptG{level,1});
            for i = 1:theta
                d1 = sum(bptG{level,i},2);
                d1(d1 == 0) = 1; % for isolated nodes
                d1_inv = d1.^(-0.5);
                D1_inv = spdiags(d1_inv, 0, N(level), N(level));
                An = D1_inv*bptG{level,i}*D1_inv;
                An = 0.5*(An + An');
                Ln_bpt{level,i} = speye(N(level)) - An;
            end
        end
    case 'asym'
        for level = 1:max_level
            N(level) = length(bptG{level,1});
            for i = 1:theta
                d1 = sum(bptG{level,i},2);
                d1(d1 == 0) = 1; % for isolated nodes
                d1_inv = d1.^(-1);
                D1_inv = spdiags(d1_inv, 0, N(level), N(level));
                An = D1_inv*(0.5*(bptG{level,i} + bptG{level,i}'));
                Ln_bpt{level,i} = speye(N(level)) - An;
            end
        end
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
c_d{1}=sgwt_cheb
y_coeff(h0,filterlen_lo,filterlen_lo+1,arange);
c_d{2}=sgwt_cheby_coeff(h1,filterlen_hi,filterlen_hi+1,arange);
c_r{1}=sgwt_cheby_coeff(g0,filterlen_hi,filterlen_hi+1,arange);
c_r{2}=sgwt_cheby_coeff(g1,filterlen_lo,filterlen_lo+1,arange);

%% Compute filter normalizations
switch GC_module
    case 'on'
        dfactor = cell(max_level,1);
        for level = 1:max_level
            for i = 1:theta
                d1 = sum(bptG{level,i},2);
                d1(d1 == 0) = 1; % for isolated nodes
                d1_sqrt = d1.^(0.5);
                dfactor{level}(:,i) = full(d1_sqrt);
            end
        end
        GC = compute_GC_Bior(Ln_bpt,dfactor, hi_d,lo_d, norm_type);
    case 'off'
        GC = cell(max_level,theta);
        for level = 1:max_level
            for j = 1:theta
                GC{level,j} = ones(N(level),2);
            end
        end
    otherwise
        disp('unknown option for GC\_module')
end

% p_lo= conv(lo_d,lo_d);
% p_hi = conv(hi_d,hi_d);
% p0 = @(x)(polyval(p_lo,x));
% p1 = @(x)(polyval(p_hi,x));
% c_p{1} = sgwt_cheby_coeff(p0,2*filterlen_lo,2*filterlen_lo+1,arange);
% c_p{2} = sgwt_cheby_coeff(p1,2*filterlen_lo,2*filterlen_lo+1,arange);
%
% for level = 1:max_level
%     for i = 1:theta
%         tk0 = zeros(N(level),1);
%         tk1 = zeros(N(level),1);
%         qk = zeros(N(level),1);
%         for iter_d = 1:30
%             v = rand(N(level),1);
%             v = double(v>0.5);
%             v = 2*v-1;
%             tk0 = tk0 + sgwt_cheby_op(v,Ln_bpt{level,i},c_p{1},arange).*v;
%             tk1 = tk1 + sgwt_cheby_op(v,Ln_bpt{level,i},c_p{2},arange).*v;
%             qk = qk + v.*v;
%         end
%         norm0_hat = mean(tk0./qk);
%         norm0(level,i) = norm0_hat^(-0.5);
%         norm1_hat = mean(tk1./qk);
%         norm1(level,i) = norm1_hat^(-0.5);
%     end
% end
% P0(P0 == 0) = 1; % for isolated nodes
% P0 = P0.^(-1);
% P0 = spdiags(P0, 0, N(1), N(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Filterbank Output at each channel
disp('Computing wavelet transform coefficients ...')
f_w = cell(max_level,1);
function_norm = zeros(max_level,1);
Channel_Name = cell(max_level,Fmax);
for level = 1:max_level
    f_w{level} = zeros(N(level)/(2^(level-1)),Fmax);
    for i = 1:Fmax
        if level == 1
            function_norm(level) = norm(f);
            tempf_w = f/function_norm(level);
        else
            function_norm(level) = norm(f_w{level-1}(Colorednodes{level-1,1},1));
            tempf_w = f_w{level-1}(Colorednodes{level-1,1},1)/function_norm(level);
        end
        for j = 1: theta
            if beta_dist{level}(i,j) == 1
                tempf_w = sgwt_cheby_op(tempf_w,Ln_bpt{level,j},c_d{1},arange);
                tempf_w = tempf_w./GC{level,j}(:,1);
                Channel_Name{level,i} = strcat(Channel_Name{level,i},'L');
            else
                tempf_w = sgwt_cheby_op(tempf_w, Ln_bpt{level,j},c_d{2},arange);
                tempf_w = tempf_w./GC{level,j}(:,2);
                Channel_Name{level,i} = strcat(Channel_Name{level,i},'H');
            end
        end
        f_w{level}(Colorednodes{level,i},i) = tempf_w(Colorednodes{level,i});
    end
end


% Plot Wavelet Coeffficents
Data_w = zeros(s_im);
x0 = 0;
y0 = 0;
dim = sqrt(N(max_level))/2;
temp_coeff = zeros(dim*dim,4);
temp_coeff(:,1) = f_w{max_level}(Colorednodes{max_level,1},1);
temp_coeff(:,2) = f_w{max_level}(Colorednodes{max_level,2},2);
temp_coeff(:,3) = f_w{max_level}(Colorednodes{max_level,3},3);
temp_coeff(:,4) = f_w{max_level}(Colorednodes{max_level,4},4);
xind = x0 + (1:dim);
yind = y0 + (1:dim);
% temp_coeff(:,2:4) = abs(temp_coeff(:,2:4));
% temp_coeff(:,1)= temp_coeff(:,1) /norm(temp_coeff(:,1));
Data_w(xind,yind) = reshape(temp_coeff(:,1),dim,dim);
y0 = y0+dim;
xind = x0 + (1:dim);
yind = y0 + (1:dim);
% temp_coeff(:,2)= temp_coeff(:,2) /norm(temp_coeff(:,2));
Data_w(xind,yind) = reshape(temp_coeff(:,2),dim,dim);
x0 = x0+dim;
y0 = 0;
xind = x0 + (1:dim);
yind = y0 + (1:dim);
% temp_coeff(:,3)= temp_coeff(:,3) /norm(temp_coeff(:,3));
Data_w(xind,yind) = reshape(temp_coeff(:,3),dim,dim);
y0 = y0+dim;
xind = x0 + (1:dim);
yind = y0 + (1:dim);
% temp_coeff(:,4)= temp_coeff(:,4) /norm(temp_coeff(:,4));
Data_w(xind,yind) = reshape(temp_coeff(:,4),dim,dim);

for level = (max_level-1):-1:1
    dim = sqrt(N(level))/2;
    temp_coeff = zeros(dim*dim,3);
    temp_coeff(:,1) = f_w{level}(Colorednodes{level,2},2);
    temp_coeff(:,2) = f_w{level}(Colorednodes{level,3},3);
    temp_coeff(:,3) = f_w{level}(Colorednodes{level,4},4);
    %     temp_coeff = abs(temp_coeff);
    x0 = 0;
    y0 = x0 + dim;
    xind = x0 + (1:dim);
    yind = y0 + (1:dim);
    %     temp_coeff(:,1)= temp_coeff(:,1) /norm(temp_coeff(:,1));
    Data_w(xind,yind) = reshape(temp_coeff(:,1),dim,dim);
    x0 = x0 +dim;
    y0 = 0;
    xind = x0 + (1:dim);
    yind = y0 + (1:dim);
    %     temp_coeff(:,2)= temp_coeff(:,2) /norm(temp_coeff(:,2));
    Data_w(xind,yind) = reshape(temp_coeff(:,2),dim,dim);
    y0 = y0 + dim;
    xind = x0 + (1:dim);
    yind = y0 + (1:dim);
    %     temp_coeff(:,3)= temp_coeff(:,3) /norm(temp_coeff(:,3));
    Data_w(xind,yind) = reshape(temp_coeff(:,3),dim,dim);
end
dim = sqrt(N(max_level))/2;
tempw = Data_w;
wav_coeffs = Data_w; % the wavelet coefficients are stored in the image format
tempw(1:dim,1:dim) = 0;
Data_w = Data_w - tempw;
% Data_w = rescale(Data_w) + rescale(abs(tempw));
Data_w = rescale(abs(tempw));
% Data_w = Data_w.^2;
scrsz = get(0,'ScreenSize');
height = scrsz(4)/1.5;
width =  scrsz(3)/2.2;
xinit = 30;
yinit = 30;
figure,
set(gcf,'Position',[xinit,yinit,width,height]);
%%
% imshow(uint8(255*Data_w))
imagesc(Data_w);
title('Wavelet Coefficients')
colormap(gray);
colorbar

%% Section 3: Non-linear Approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_w = find_kbest(f_w, nnz_factor, max_level, Colorednodes); % threshold for nCoeffs best wavelet coefficients

%% Section 4: Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct Signals in Original Graph Domain by using Coefficients
% Compute Filterbank Output at each channel
disp('Reconstructed signals using single channel coefficients ...')
f_hat1 = zeros(N(1),1);
for d_level = max_level
    f_hat = f_w;
    for level = d_level: -1 :1
        % Channel_Name = cell(Fmax,1);
        for i = 1:Fmax
            tempf_hat = f_hat{level}(:,i);
            for j = theta: -1: 1
                if beta_dist{level}(i,j) == 1
                    tempf_hat = GC{level,j}(:,1).*tempf_hat;
                    tempf_hat = sgwt_cheby_op(tempf_hat,Ln_bpt{level,j},c_r{1},arange);
                else
                    tempf_hat = GC{level,j}(:,2).*tempf_hat;
                    tempf_hat = sgwt_cheby_op(tempf_hat,Ln_bpt{level,j},c_r{2},arange);
                end
            end
            f_hat{level}(:,i) = tempf_hat;
        end
        f_hat{level} = sum(f_hat{level},2)*function_norm(level);
        %         f_hat{level} = f_hat{level}(:,1);
        
        if level >1
            f_hat{level-1}(Colorednodes{level-1,1},1) = f_hat{level};
        end
    end
    f_hat1= f_hat{level};
end
figure,
subplot(2,1,1)
c_image = reshape(f,s_im,s_im); % original image
imagesc(c_image)
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('Original Image')
colorbar
subplot(2,1,2)
c_image = reshape(f_hat1,s_im,s_im); % approximation coeffs
imagesc(c_image)
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
S = sprintf('%s%1.2f%s',' reconstruction from ', nnz_factor,' highpass coefficients');
title(S)
colormap(gray)
colorbar

% % Plot reconstructed signals
% rowsize = 3;
% colsize = ceil(max_level/rowsize)+1;
% gap = 0.1;
% dr = (1 - (rowsize+1)*gap)/rowsize;
% dc = (1 - (colsize+1)*gap)/colsize;
% s_im = sqrt(N(1));
% xinit = width+30;
% yinit = 30;
% figure,
% set(gcf,'Position',[xinit,yinit,width,height]);
% subplot('Position',[gap (1-gap-dc) dr dc])
% c_image = reshape(f,s_im,s_im); % original image
% imagesc(c_image)
% set(gca,'Xtick',[]);
% set(gca,'Ytick',[]);
% title('Original Image')
% for level = 1:max_level
%     col = floor((level -1)/rowsize) + 1;
%     row = level - (col -1)*rowsize;
%     yinit =   1 - (gap+dc)*(col+1);
%     xinit = gap + (gap + dr)*(row-1);
%     subplot('Position',[xinit yinit dr dc])
%     c_image = reshape(f_hat1(:,level),s_im,s_im); % approximation coeffs
%     imagesc(c_image)
%     set(gca,'Xtick',[]);
%     set(gca,'Ytick',[]);
%     S1 = repmat('L',1,level);
%     S = sprintf('%s%s',' reconstruction from ', S1);
%     title(S)
%     colormap(gray)
% end

%% Section 5: Summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uLam = 0:1/N:2;

h0 = sgwt_cheby_eval(uLam,c_d{1},arange);
h1 = sgwt_cheby_eval(uLam,c_d{2},arange);
g0 = sgwt_cheby_eval(uLam,c_r{1},arange);
g1 = sgwt_cheby_eval(uLam,c_r{2},arange);
figure,
plot(uLam,h0);
hold on
plot(uLam,h1,'r');
title('Spectral Kernels-Analysis');
legend('Low-pass', 'HighPass');
figure,
plot(uLam,g0);
hold on
plot(uLam,g1,'r');
title('Spectral Kernels-Synthesis');
legend('Low-pass', 'HighPass');
% MSE_kernel = mean((g1-g(uLam)).^2);
Err = 1 - 0.5*(h1.*g1 + h0.*g0);
Err = max(abs(Err));
S = sprintf('%s%f','The max Error in kernel-approximation = ' ,Err);
disp(S);
% f1 = repmat(f,1,max_level);
MSE = (f - f_hat1).^2;
MSE = sum(MSE);
SNR = 10*log10(norm(f)^2./MSE);
S = sprintf('%s%f','The SNR of reconstruction is = ' ,SNR);
disp(S);
end
%% Function for non-linear approximation
function f = find_kbest(f_w, nnz, d_level,Colorednodes)
if nnz < 1
    lin_vec = [];
    for level = 1:d_level
        lin_vec = [lin_vec; f_w{level}(Colorednodes{level,2},2); f_w{level}(Colorednodes{level,3},3); f_w{level}(Colorednodes{level,4},4)];
    end
    nCoeffs = floor(length(lin_vec)*nnz);
    lin_vec = sort(abs(lin_vec(:)),'descend');
    thresh = lin_vec(nCoeffs+1);
    for level = 1:d_level
        temp = f_w{level}(:,2:4);
        temp(abs(temp) <thresh) = 0;
        %     temp(isolated{level},:) = f_w{level}(isolated{level},2:4);
        f_w{level}(:,2:4) = temp;
    end
end
f = f_w;
end