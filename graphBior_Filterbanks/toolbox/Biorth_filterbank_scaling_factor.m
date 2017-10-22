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
%   Date last modified:	02/01/2012 Sunil K. Narang
%
%   Description:
%   This file generates QMF filter-banks on the 8-connected graph formulation
%   of 2D digital images as proposed in the paper:
%% "S. K. Narang and Antonio Ortega, "Perfect Reconstruction Two-Channel
%%  Wavelet Filter-Banks For Graph Structured Data",
%  IEEE TSP also avaliable as Tech. Rep. arXiv:1106.3693v3, Dec 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bio_Scale] = Biorth_filterbank_scaling_factor(varargin)
len = size(varargin,2);
switch len
    case 0,
        filename = 'lena.png';
        filetype = 'png';
        opt = struct('max_level',4,'vanishing_moments',[5,7],'nnz_factor',1,'edgemap', 0);
    case 2,
        filename = varargin{1};
        filetype = varargin{2};
        opt = struct('max_level',2,'vanishing_moments',[7,9],'nnz_factor',1,'edgemap', 0);
    case 3,
        filename = varargin{1};
        filetype = varargin{2};
        opt = varargin{3};
    otherwise,
        disp('Invalid # of input arguments. Enter inputs as (filename,filetype,option)');
        return;
end

addpath Graph_Generators/
addpath Graph_kernels/
addpath sgwt_toolbox/
addpath toolbox/

% Parameters
max_level = opt.max_level; % number of decomposition levels
if isempty(max_level)
    max_level =2;
end
vanishing_moments = opt.vanishing_moments; % filter length
if isempty(vanishing_moments)
    filterlen = [7,9];
end
nnz_factor = opt.nnz_factor; % fraction of non-zero coefficient
if isempty(nnz_factor)
    nnz_factor = 1;%(4^(max_level-2)); 
end
if isempty(filename)
     filename = 'Lichtenstein.png';
     filetype = 'png';
end 
edgemap =opt.edgemap; % uses edge-map if 1, use regular 8-connected graph otherwise
if isempty(edgemap)
    edgemap = 0;%(4^(max_level-2)); 
end
theta = 2; % number of bipartite graphs
Fmax = 2^theta; % Graph Coloring
N = zeros(max_level,1); %size of bipartite graphs at each level
norm_type = 'asym'; % use 'sym'for symmetric Normalized Laplacian matrix
S = sprintf('%s%d%s%d%s%d%s%f%s%d', '# decomposition levels = ', max_level,'.  Biorthtype = ', vanishing_moments(1),' / ',vanishing_moments(2), '.   nnz_factor = ',nnz_factor,'.   Using edgemap = ', edgemap);
disp(S);


%% Section 1: Image Graph Formulation
% Graph Signal
Data = imread(filename,filetype);
if length(size(Data)) == 3
    Data = rgb2gray(Data);
end
Data = double(Data);%/255;
[m n] = size(Data);
m = min([m,n]);
s_im = floor(m/2^max_level)*2^max_level;
% s_im = 10;
Data = Data(1:s_im,1:s_im);
f_or=Data(:);
f_or = f_or/norm(f_or);
[bptG, Colorednodes, beta_dist, loc] = image_graphs_multi(Data,max_level,edgemap);  % image graphs

%noise_signal=Noise(:);
%noise_signal=noise_signal/norm(f);

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
% design a low-pass kernel
Nd = vanishing_moments(1);
Nr = vanishing_moments(2);
% S = sprintf('%s%d%s ' , ' Computing a ',filterlen,' ^th order approximation of Meyer kernel');
% disp(S)
[hi_d,lo_d] = biorth_kernel(Nd,Nr); % Nd zeros for lowpass Nr zeros for highpass
filterlen_lo = length(roots(lo_d));
filterlen_hi = length(roots(hi_d));
h0 = @(x)(polyval(lo_d,x));
h1 = @(x)(polyval(hi_d,x));
g0 = @(x)(polyval(hi_d,2 - x));
g1 = @(x)(polyval(lo_d,2 - x));
arange = [0 2];
c_d{1}=sgwt_cheby_coeff(h0,filterlen_lo,filterlen_lo+1,arange);
c_d{2}=sgwt_cheby_coeff(h1,filterlen_hi,filterlen_hi+1,arange);
c_r{1}=sgwt_cheby_coeff(g0,filterlen_hi,filterlen_hi+1,arange);
c_r{2}=sgwt_cheby_coeff(g1,filterlen_lo,filterlen_lo+1,arange);

%%% Compute filter normalizations
p_lo= conv(lo_d,lo_d);
p_hi = conv(hi_d,hi_d);
p0 = @(x)(polyval(p_lo,x));
p1 = @(x)(polyval(p_hi,x));
c_p{1} = sgwt_cheby_coeff(p0,2*filterlen_lo,2*filterlen_lo+1,arange);
c_p{2} = sgwt_cheby_coeff(p1,2*filterlen_lo,2*filterlen_lo+1,arange);
for level = 1:max_level
    for i = 1:theta
%         tempf_w = speye(N(level));
%         tempP0= sgwt_cheby_op(tempf_w,Ln_bpt{level,i},c_p{i},arange);
%         P0{level,i} = diag(tempP0).^(0.5);
        P0{level,i} = diag(speye(N(level)));
    end
end

max_iter = 20;
Bio_Scale = zeros(max_level,Fmax);
for iter = 1:max_iter
    Noise = randn(s_im,s_im);
    %Data=Data + Noise;
    Data = Noise;
    f = Data(:);
    disp('Computing wavelet transform coefficients ...')
    f_w = cell(max_level,1);
    Channel_Name = cell(max_level,Fmax);
    for level = 1:max_level
        f_w{level} = zeros(N(level)/(2^(level-1)),Fmax);
        for i = 1:Fmax
            if level == 1
                tempf_w = f;
            else
                tempf_w = f_w{level-1}(Colorednodes{level-1,1},1);
            end
            for j = 1: theta
                if beta_dist{level}(i,j) == 1
                    tempf_w = sgwt_cheby_op(tempf_w,Ln_bpt{level,j},c_d{1},arange);
                    tempf_w = (P0{level,j}.^(-1)).*tempf_w;
                    Channel_Name{level,i} = strcat(Channel_Name{level,i},'L');
                else
                    tempf_w = sgwt_cheby_op(tempf_w,Ln_bpt{level,j},c_d{2},arange);
                    tempf_w = (P0{level,j}.^(-1)).*tempf_w;
                    Channel_Name{level,i} = strcat(Channel_Name{level,i},'H');
                end
            end
            f_w{level}(Colorednodes{level,i},i) = tempf_w(Colorednodes{level,i});
            Bio_Scale(level,i) = Bio_Scale(level,i) + var(tempf_w(Colorednodes{level,i}));
        end
    end
end
Bio_Scale = Bio_Scale/max_iter;
Bio_Scale = Bio_Scale.^(0.5);
% %% Section 3: Non-linear Approximation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % f_w = find_kbest(f_w, nnz_factor, max_level, Colorednodes); % threshold for nCoeffs best wavelet coefficients
% 
% Noise = randn(s_im,s_im)*20;
% %Data=Data + Noise;
% Data = Noise;
% f = Data(:);
% disp('Computing wavelet transform coefficients ...')
% f_w = cell(max_level,1);
% Channel_Name = cell(max_level,Fmax);
% for level = 1:max_level
%     f_w{level} = zeros(N(level)/(2^(level-1)),Fmax);
%     for i = 1:Fmax
%         if level == 1
%             tempf_w = f;
%         else
%             tempf_w = f_w{level-1}(Colorednodes{level-1,1},1);
%         end
%         for j = 1: theta
%             if beta_dist{level}(i,j) == 1
%                 tempf_w = sgwt_cheby_op(tempf_w,Ln_bpt{level,j},c_d{1},arange);
%                 Channel_Name{level,i} = strcat(Channel_Name{level,i},'L');
%             else
%                 tempf_w = sgwt_cheby_op(tempf_w,Ln_bpt{level,j},c_d{2},arange);
%                 Channel_Name{level,i} = strcat(Channel_Name{level,i},'H');
%             end
%         end
%         f_w{level}(Colorednodes{level,i},i) = tempf_w(Colorednodes{level,i});
%     end
% end
% % Plot Wavelet Coeffficents
% Data_w = zeros(s_im);
% x0 = 0;
% y0 = 0;
% dim = sqrt(N(max_level))/2;
% temp_coeff = zeros(dim*dim,4);
% temp_coeff(:,1) = f_w{max_level}(Colorednodes{max_level,1},1);
% temp_coeff(:,2) = f_w{max_level}(Colorednodes{max_level,2},2);
% temp_coeff(:,3) = f_w{max_level}(Colorednodes{max_level,3},3);
% temp_coeff(:,4) = f_w{max_level}(Colorednodes{max_level,4},4);
% 
% xind = x0 + (1:dim);
% yind = y0 + (1:dim);
% % temp_coeff(:,2:4) = abs(temp_coeff(:,2:4));
% % temp_coeff(:,1)= temp_coeff(:,1) /norm(temp_coeff(:,1));
% %Data_w(xind,yind) = reshape(temp_coeff(:,1),dim,dim)./(Bio_Scale(1)*Bio_Scale(1)^(max_level-1));
% Data_w(xind,yind) = reshape(temp_coeff(:,1),dim,dim)/Bio_Scale(max_level,1);
% y0 = y0+dim;
% xind = x0 + (1:dim);
% yind = y0 + (1:dim);
% % temp_coeff(:,2)= temp_coeff(:,2) /norm(temp_coeff(:,2));
% %Data_w(xind,yind) = reshape(temp_coeff(:,2),dim,dim)./(Bio_Scale(2)*Bio_Scale(1)^(max_level-1));
% Data_w(xind,yind) = reshape(temp_coeff(:,2),dim,dim)/Bio_Scale(max_level,2);
% x0 = x0+dim;
% y0 = 0;
% xind = x0 + (1:dim);
% yind = y0 + (1:dim);
% % temp_coeff(:,3)= temp_coeff(:,3) /norm(temp_coeff(:,3));
% %Data_w(xind,yind) = reshape(temp_coeff(:,3),dim,dim)./(Bio_Scale(3)*Bio_Scale(1)^(max_level-1));
% Data_w(xind,yind) = reshape(temp_coeff(:,3),dim,dim)/Bio_Scale(max_level,3);
% y0 = y0+dim;
% xind = x0 + (1:dim);
% yind = y0 + (1:dim);
% % temp_coeff(:,4)= temp_coeff(:,4) /norm(temp_coeff(:,4));
% %Data_w(xind,yind) = reshape(temp_coeff(:,4),dim,dim)./(Bio_Scale(4)*Bio_Scale(1)^(max_level-1));
% Data_w(xind,yind) = reshape(temp_coeff(:,4),dim,dim)/Bio_Scale(max_level,4);
% 
% for level = (max_level-1):-1:1
%     dim = sqrt(N(level))/2;
%     temp_coeff = zeros(dim*dim,3);
%     temp_coeff(:,1) = f_w{level}(Colorednodes{level,2},2);
%     temp_coeff(:,2) = f_w{level}(Colorednodes{level,3},3);
%     temp_coeff(:,3) = f_w{level}(Colorednodes{level,4},4);
%     
%     
%     %     temp_coeff = abs(temp_coeff);
%     x0 = 0;
%     y0 = x0 + dim;
%     xind = x0 + (1:dim);
%     yind = y0 + (1:dim);
%     %     temp_coeff(:,1)= temp_coeff(:,1) /norm(temp_coeff(:,1));
%     Data_w(xind,yind) = reshape(temp_coeff(:,1),dim,dim)/Bio_Scale(level,2);
%     x0 = x0 +dim;
%     y0 = 0;
%     xind = x0 + (1:dim);
%     yind = y0 + (1:dim);
%     %     temp_coeff(:,2)= temp_coeff(:,2) /norm(temp_coeff(:,2));
%     Data_w(xind,yind) = reshape(temp_coeff(:,2),dim,dim)/Bio_Scale(level,3);
%     y0 = y0 + dim;
%     xind = x0 + (1:dim);
%     yind = y0 + (1:dim);
%     %     temp_coeff(:,3)= temp_coeff(:,3) /norm(temp_coeff(:,3));
%     Data_w(xind,yind) = reshape(temp_coeff(:,3),dim,dim)/Bio_Scale(level,4);
% end
% dim = sqrt(N(max_level))/2;
% tempw = Data_w;
% wav_coeffs = Data_w; % the wavelet coefficients are stored in the image format
% %figure,
% %imagesc(wav_coeffs)
% %colormap(gray);
% % tempw(1:dim,1:dim) = 0;
% % Data_w = Data_w - tempw;
% % Data_w = rescale(Data_w) + rescale(abs(100*tempw));
% Data_w = rescale(Data_w);
% scrsz = get(0,'ScreenSize');
% height = scrsz(4)/1.5;
% width =  scrsz(3)/2.2;
% xinit = 30;
% yinit = 30;
% figure,
% set(gcf,'Position',[xinit,yinit,width,height]);
% imshow(uint8(255*Data_w))
% %imagesc(Data_w)
% title('Wavelet Coefficients Before Denoising')
% colormap(gray);
