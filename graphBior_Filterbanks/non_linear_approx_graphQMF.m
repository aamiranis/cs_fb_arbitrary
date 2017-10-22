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
function Output = non_linear_approx_graphQMF(varargin)
clc
% close all
disp('This function computes wavelet coefficients of any Image using a 2-dim separable graph-QMF filterbank.')
disp('See help for input options');
addpath Graph_Generators/
addpath Graph_kernels/
addpath sgwt_toolbox/
addpath toolbox/
addpath exportfig/
addpath ImageQualityMeasures/

p = inputParser;
p.addOptional('filename','ballet.png',@(f) exist(f,'file'));
p.addParamValue('maxLevel',3,@isscalar);
p.addParamValue('filterLen',18,@isscalar);
p.addParamValue('LType','asym',@isstr);
p.addParamValue('edgemap',0,@islogical);
p.parse(varargin{:});
opt = p.Results;

% Set parameters
filename = opt.filename;
tempI = imfinfo(filename);
filetype = tempI.Format;
max_level = opt.maxLevel; % number of decomposition levels
filterlen = opt.filterLen; % filter length
norm_type = opt.LType;
edgemap = opt.edgemap;
S = sprintf('%s%d\n%s%d\n%s%d\n%s%s', '# decomposition levels = ', max_level,'Filterlength = ', filterlen, 'Using edgemap = ', edgemap,'Type of Laplacian = ',norm_type);
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
% s_im = 512;
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
% design QMF kernels
S = sprintf('%s%d%s ' , ' Computing a ',filterlen,' ^th order approximation of Meyer kernel');
disp(S)
%g = @(x)(ideal_kernel(x)); % Ideal Wavelet Kernel
g_low = @(x)(meyer_kernel(x));  % Meyer Wavelet Kernel
g_high =@(x)(meyer_kernel(2-x));
%g = @(x)(johnston_kernel(x,filterlen));  % Meyer Wavelet Kernel
arange = [0 2];
c_d{1}=sgwt_cheby_coeff(g_low,filterlen,filterlen+1,arange);
c_d{2}=sgwt_cheby_coeff(g_high,filterlen,filterlen+1,arange);
c_r{1} = c_d{1};
c_r{2} = c_d{2};

GC_module = 'off';
%% Compute filter normalizations
disp('Computing Gain in each channel. This may take some time ...');
switch GC_module
    case 'on'
        switch norm_type
            case 'sym'
                GC = compute_GC_Bior(Ln_bpt,[], hi_d,lo_d, norm_type);
                %                 GC_exact = compute_GC_Bior_exact(Ln_bpt,[], hi_d,lo_d, norm_type);
            case 'asym'
                dfactor = cell(max_level,1);
                Ln_sym = cell(max_level,theta);
                for level = 1:max_level
                    for i = 1:theta
                        A2 = double(bptG{level,i}>0);
                        d1 = sum(bptG{level,i},2);
                        d1(d1 == 0) = 1; % for isolated nodes
                        d1_inv = d1.^(-1);
                        numnbr = sum(A2,2);
                        numnbr(numnbr == 0) = 1; % for isolated nodes
                        numnbr_inv = numnbr.^(-1);
                        dfactor{level}(:,i) = (d1_inv.*numnbr_inv).*smvp(bptG{level,i},d1);
                        d1_inv = d1.^(-0.5);
                        D1_inv = spdiags(d1_inv, 0, N(level), N(level));
                        An = D1_inv*bptG{level,i}*D1_inv;
                        An = 0.5*(An + An');
                        Ln_sym{level,i} = speye(N(level)) - An;
                    end
                end
                GC = compute_GC_Bior(Ln_sym,dfactor, hi_d,lo_d, norm_type);
        end
    case 'off'
        %         switch norm_type
        %             case 'sym'
        %                 GC = compute_GC_Bior(Ln_bpt,[], hi_d,lo_d, norm_type);
        %                 %                 GC_exact = compute_GC_Bior_exact(Ln_bpt,[], hi_d,lo_d, norm_type);
        %             case 'asym'
        %                 Ln_sym = cell(max_level,theta);
        %                 for level = 1:max_level
        %                     for i = 1:theta
        %                         d1 = sum(bptG{level,i},2);
        %                         d1(d1 == 0) = 1; % for isolated nodes
        %                         d1_inv = d1.^(-0.5);
        %                         D1_inv = spdiags(d1_inv, 0, N(level), N(level));
        %                         An = D1_inv*bptG{level,i}*D1_inv;
        %                         An = 0.5*(An + An');
        %                         Ln_sym{level,i} = speye(N(level)) - An;
        %                     end
        %                 end
        %                 GC = compute_GC_Bior(Ln_sym,[], hi_d,lo_d, 'sym');
        %         end
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
function_norm = ones(max_level,1);
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
Output{1} = Data_w;
tempw = Data_w(1:dim,(dim+1):s_im);
tempw1 = Data_w((dim+1):s_im,:);
tempw = [tempw(:); tempw1(:)];
% hist(abs(tempw),64);
% set(gcf, 'Color', 'w');
% V = axis;
% V(2) = 9*10^-3;
% axis(V);
% outpfile = sprintf('%s%s%s%s%d%s%d%s%d%s%s%s%s%s','Results/',filename,'_hist','_FL',filterlen,'_edgemap',edgemap,'_',norm_type,'_GC',GC_module,'.png');
% export_fig(outpfile,'-png','-m2');
norm0 = length(find(abs(tempw) >10^-6));
Output{2} = [norm0 norm(tempw,1) norm(tempw,2) norm(tempw,Inf) norm(tempw,-Inf)];
Output{3} = abs(tempw);
% tempw = Data_w;
% tempw(1:dim,1:dim) = 0;
% Data_w = Data_w - tempw;
% Data_w = rescale(Data_w) + rescale(abs(tempw));
% % Data_w = Data_w.^2;
% outpfile = sprintf('%s%s%s%s%d%s%d%s%s%s%s%s','Results/',filename,'_wavelet_coeffs','_FL',filterlen,'_levels',max_level,'_',norm_type,'_GC',GC_module,'.png');
% imwrite(Data_w,outpfile,'png');


origImg = reshape(f,s_im,s_im);
nnz = 0:0.01:0.18;
MSE = zeros(1,length(nnz));
SSIM = zeros(1,length(nnz));
PSNR = zeros(1,length(nnz));
SNR = zeros(1,length(nnz));
f_w1 = f_w;
for iter = 1:length(nnz)
    nnz_factor = nnz(iter);
    %         f_w = find_kbest(f_w1, nnz_factor, max_level, Colorednodes); % threshold for nCoeffs best wavelet coefficients
    %     wav_coeffs = wavelet_ordering(f_w,s_im, Colorednodes);
    %     figure,imagesc(wav_coeffs);
    %     colormap(gray);
    %     title('without scaling')
    
    [f_w thresh(iter)] = find_kbest_biorth(f_w1, nnz_factor, max_level, Colorednodes); % threshold for nCoeffs best wavelet coefficients
    %     for level = 1:max_level
    %         f_w{level} = f_w1{level}-f_w{level};
    %     end
    %     f_w2 = f_w;
    %     for level = 1:max_level
    %         for F = 1:Fmax
    %             f_w2{level}(:,F) = f_w2{level}(:,F)/Bio_Scale(level,F);
    %         end
    %     end
    %     wav_coeffs = wavelet_ordering(f_w2,s_im, Colorednodes);
    %     if iter == 1
    %         cLim = [min(wav_coeffs(:)), max(wav_coeffs(:))];
    %     end
    %     figure,imagesc(wav_coeffs);
    %     colormap(gray);
    %     colorbar
    %     title('with scaling')
    
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
    
    %% Section 5: Summary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    distImg = reshape(f_hat1,s_im,s_im);
    if iter == 2
        %         outpfile = sprintf('%s%s%s%s%d%s%d%s%d%s%s%s%s%s','Results/',filename,'_outp','_FL',filterlen,'_levels',max_level,'_NNZ',nnz_factor*100,'%_',norm_type,'_GC',GC_module,'.png');
        imagesc(distImg);
        colorbar
        colormap(gray)
        set(gca,'xTick',[])
        set(gca,'yTick',[])
    end
    %     figure,imagesc(distImg);
    %     colormap(gray)
    %     colorbar
    %     Mean Square Error
    MSE(iter) = MeanSquareError(origImg, distImg);
    % disp('Mean Square Error = ');
    % disp(MSE);
    
    %Peak Signal to Noise Ratio
    PSNR(iter) = PeakSignaltoNoiseRatio(origImg, distImg);
    % disp('Peak Signal to Noise Ratio = ');
    % disp(PSNR);
    
    SSIM(iter) = ssim(origImg,distImg);
    
    SE = (f - f_hat1).^2;
    SE = sum(SE);
    if(SE > 0)
        SNR(iter) = 10*log(norm(f)^2./SE)/log(10);
    else
        SNR(iter) = 99;
    end
end
Output{4} = struct('SNR',SNR,'PSNR',PSNR, 'SSIM',SSIM,'MSE', MSE);