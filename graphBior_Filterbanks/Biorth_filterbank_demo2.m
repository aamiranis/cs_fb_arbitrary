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
%   This file generates QMF filter-banks on the Minnesota traffic graph
%   example given in the paper:
%% "S. K. Narang and Antonio Ortega, "Perfect Reconstruction Two-Channel
%%  Wavelet Filter-Banks For Graph Structured Data",
%  IEEE TSP also avaliable as Tech. Rep. arXiv:1106.3693v3, Dec 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [wav_coeffs, channel_info] = Biorth_filterbank_demo2()
clear all
clc
close all
disp('***********Welcome to Demo 2 of Graph-Bior*********')
disp('In this demo we implement a 2-dim separable graph-Bior filterbank on Minnesota traffic graph')
profile on
addpath Graph_Generators/
addpath Graph_kernels/
addpath sgwt_toolbox/
addpath toolbox/
addpath Bpt_Decomp_Algos/
addpath sgwt_toolbox/
addpath Datasets/

%% Parameters
max_level = 1; % number of decomposition levels
nnz_factor = 0.04;
filterlen = 10; % filterlength
norm_type = 'asym';
GC_module = 'off';
S = sprintf('%s%d\n%s%d\n%s%s\n%s%f\n%s%s', '# decomposition levels = ', max_level,'Filterlength = ', filterlen,'GC module is ',GC_module, 'Fraction of high-pass coefficients used = ',nnz_factor,'Type of Laplacian = ',norm_type);
disp(S);
%% Section 1: Graph Formulation

% Minnesota traffic graph
[A xy] = Minnesota_traffic_graph(); % Minnesota Traffic Graph
% [A xy] = random_bipartite_graph(); % random bipartite graphs G(n,m,p)
% N = 500;
% [A xy] = generate_graph('sensor',N,0);
% load 685_bus
% A = Problem.A;
% load bus_ex_coord.x
% A = double(abs(A)>0);
% xy = bus_ex_coord;
%% Section 2: Bipartite subgraph decomposition

% Graph Coloring: F is the output coloring vector
%%%%%%%%%%%%%
% F = BSC(A); % Backtracking Sequential Coloring Algorithm: Exact but slow
% F = LF(A);   % Greedy Largest Degree First Coloring: Fast but inaccurate
% F = DSATUR(A); % Greedy Coloring based on improved heuristic: Moderate speed and accuracy
load min_coloring %  Preloaded Coloring for Minnesota graph
%%%%%%%%%%%%%%

% generate the downsampling functions
N = length(A);
[beta bptG beta_dist Colorednodes]= harary_decomp1(A,F);
theta = size(beta,2);
Fmax = size(beta_dist,1);

% f = mean(xy.^2,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 4: Filterbank implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Normalized Laplacian Matrices for Each Bpt graph
disp('Computing normalized Laplacian matrices for each subgraph...');
Ln_bpt = cell(1,theta);
switch norm_type
    case 'sym'
        for i = 1:theta
            %     bptG(:,:,i) = 0.5*(bptG(:,:,i) + bptG(:,:,i)');
            d1 = sum(bptG(:,:,i),2);
            d1(d1 == 0) = 1; % for isolated nodes
            d1_inv = d1.^(-0.5);
            D1_inv = diag(d1_inv);
            An = D1_inv*bptG(:,:,i)*D1_inv;
            An = 0.5*(An + An');
            %     d1_inv = d1.^(-1);
            %     D1_inv = diag(d1_inv);
            %     An = D1_inv*bptG(:,:,i);
            Ln_bpt{1,i} = eye(N) - An;
        end
        
    case 'asym'
        for i = 1:theta
            bptG(:,:,i) = 0.5*(bptG(:,:,i) + bptG(:,:,i)');
            d1 = sum(bptG(:,:,i),2);
            d1(d1 == 0) = 1; % for isolated nodes
            d1_inv = d1.^(-1);
            D1_inv = diag(d1_inv);
            An = D1_inv*bptG(:,:,i);
            Ln_bpt{1,i} = eye(N) - An;
        end
        
    otherwise
        disp('Unknown normalization option')
        return;
end



% design a low-pass kernel
% design biorthogonal kernels
N_lo = floor(filterlen/2) + mod(filterlen,2);
N_hi = N_lo;
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


%% Compute filter normalizations
% switch GC_module
%     case 'on'
%         disp('Computing Gain in each channel. This may take some time ...');
%         switch norm_type
%             case 'sym'
%                 GC = compute_GC_Bior(Ln_bpt,[], hi_d,lo_d, norm_type);
% %                 GC_exact = compute_GC_Bior_exact(Ln_bpt,[], hi_d,lo_d, norm_type);
%             case 'asym'
%                 dfactor = cell(max_level,1);
%                 for level = 1:max_level
%                     for i = 1:theta
%                         d1 = sum(bptG{level,i},2);
%                         d1(d1 == 0) = 1; % for isolated nodes
%                         d1_inv = d1.^(-1);
%                         dfactor{level}(:,i) = d1_inv.*smvp(bptG{level,i},d1);
%                     end
%                 end
%                 GC = compute_GC_Bior(Ln_bpt,dfactor, hi_d,lo_d, norm_type);
%         end
%     case 'off'
%         disp('Skipping GC module');
%         GC = cell(max_level,theta);
%         for level = 1:max_level
%             for j = 1:theta
%                 GC{level,j} = ones(N(level),2);
%             end
%         end
%     otherwise
%         disp('unknown option for GC\_module')
% end
%
switch GC_module
    case 'on'
        dfactor = cell(max_level,1);
        for level = 1:max_level
            for i = 1:theta
                d1 = sum(bptG(:,:,i),2);
                d1(d1 == 0) = 1; % for isolated nodes
                d1_sqrt = d1.^(0.5);
                dfactor{level}(:,i) = full(d1_sqrt);
            end
        end
        GC = compute_GC_Bior_graph(Ln_bpt,beta,dfactor, hi_d,lo_d, norm_type);
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


% S = sprintf('%s%d%s ' , ' Computing a ',filterlen,' ^th order approximation of Meyer kernel');
%% Section 3: Graph Signal
% if exist('min_graph_signal.m','file')
load min_graph_signal % Preloaded graph signal for Minnesota traffic
% else
% graph
if ~exist('f','var') % or 
    D = diag(sum(A,2));
    d = sum(A,2);
    d(d == 0) = 1;
    d_inv = d.^(-0.5);
    D_inv = diag(d_inv);
    Ln = D_inv*A*D_inv;
    Ln = 0.5*(Ln + Ln');
    [U Lam] = eigs(Ln);
    f = U(:,4);
    f = f - mean(f);        
    f(f<0) = -1;
    f(f>=0) = 1;
    f = double(f);
%     f(f==0) = -1;
    save min_graph_signal f
end
% f1 = f;
% r = randn(size(f));
% f = f+0.1*r;
% MSE_orig = norm(f1 - sum(f,2))^2;
% PSNR_orig = 10*log10(norm(f1)^2/MSE_orig);
% S = sprintf('%s%f','The SNR of noisy signal = ' ,PSNR_orig);
% disp(S);
% f = ones(N,1);
% P_l = U(:,1:N_lo)*U(:,1:N_lo)';
% f = P_l*f;
% disp(S)

disp('Computing wavelet transform coefficients ...')
function_norm = zeros(max_level,1);
Channel_Name = cell(max_level,Fmax);
f_w = zeros(N,Fmax);
for i = 1:Fmax
%         function_norm(1) = norm(f);
%         tempf_w = f/function_norm(1);
          tempf_w = f;
    for j = 1: theta
        if beta_dist(i,j) == 1
            tempf_w = sgwt_cheby_op(tempf_w,Ln_bpt{j},c_d{1},arange);
            tempf_w = tempf_w./GC{1,j}(:,1);
            Channel_Name{1,i} = strcat(Channel_Name{1,i},'L');
        else
            tempf_w = sgwt_cheby_op(tempf_w, Ln_bpt{1,j},c_d{2},arange);
            tempf_w = tempf_w./GC{1,j}(:,2);
            Channel_Name{1,i} = strcat(Channel_Name{1,i},'H');
        end
    end
    if ~isempty(Colorednodes{i})
        f_w(Colorednodes{i},i) = tempf_w(Colorednodes{i});
    end
end

wav_coeffs = sum(f_w,2);
channel_info(Fmax).name =[];
for i = 1:Fmax
    channel_info(i).name = Channel_Name{i};
    channel_info(i).nodes = Colorednodes{i};
end
channel_info(3).name = 'HL';
f_w1 = f_w;
details = [f_w(Colorednodes{2},2) ; f_w(Colorednodes{4},4)];
[~,I]= sort(abs(details),'descend');
cutoff = ceil(N*nnz_factor);
details(I(cutoff:end)) = 0;
f_w(Colorednodes{2},2) = details(1:length(Colorednodes{2}));
f_w(Colorednodes{4},4) = details((length(Colorednodes{2})+1):end);

%% Section 5: Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct Signals in Original Graph Domain by using Coefficients
% Compute Filterbank Output at each channel
disp('Reconstructed signals using single channel coefficients ...')
f_hat = f_w;
% Channel_Name = cell(Fmax,1);
for i = 1:Fmax
    tempf_hat = f_hat(:,i);
    for j = theta: -1: 1
        if beta_dist(i,j) == 1
            tempf_hat = GC{1,j}(:,1).*tempf_hat;
            tempf_hat = sgwt_cheby_op(tempf_hat,Ln_bpt{1,j},c_r{1},arange);
        else
            tempf_hat = GC{1,j}(:,2).*tempf_hat;
            tempf_hat = sgwt_cheby_op(tempf_hat,Ln_bpt{1,j},c_r{2},arange);
        end
    end
    f_hat(:,i) = tempf_hat;
end
f_hat = sum(f_hat,2);%*function_norm(1);
%         f_hat{level} = f_hat{level}(:,1);

f_hat1= f_hat;
    
%% Section 6: Summary
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uLam = 0:1/N:2;
% % if strcmp(filter_type,'exact')
% %     g1 = g(uLam);
% % else
%     h0 = sgwt_cheby_eval(uLam,c_d{1},arange);
%     h1 = sgwt_cheby_eval(uLam,c_d{2},arange);
%     g0 = sgwt_cheby_eval(uLam,c_r{1},arange);
%     g1 = sgwt_cheby_eval(uLam,c_r{2},arange);
% % end
% % MSE_kernel = mean((g1-g(uLam)).^2);
% Err = 1 - 0.5*(h1.*g1 + h0.*g0);
% Err = max(abs(Err));
% S = sprintf('%s%f','The max Error in kernel-approximation = ' ,Err);
% disp(S);
MSE= mean((f - f_hat).^2);
P_rms = mean(f.*f);
SNR = 10*log10(P_rms/MSE);
S = sprintf('%s%f','The MSE between overall reconstructed graph signal and original signal = ' ,MSE);
disp(S);
S = sprintf('%s%f','The SNR of reconstruction = ' ,SNR);
disp(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 7: Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bipartite Decomposition
if ~exist('xy','var')
    [x y] = draw_dot(A);
    xy(:,1) = x'; xy(:,2) = y';
end
scrsz = get(0,'ScreenSize');
height = scrsz(4)/5;
width = scrsz(3)/5;

figure,
gplot(A,xy,'r.-')
title('Original Graph')
set(gcf,'Position',[30,30,width,height]);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

for i = 1:theta
    figure,
    S1 = find(beta(:,i) == 1);
    S2 = find(beta(:,i) == -1);
    plot(xy(S1,1),xy(S1,2),'bs');
    hold on
    plot(xy(S2,1),xy(S2,2),'ro');
    gplot(bptG(:,:,i),xy,'r-');
    S = sprintf('%s%d','bipartite subgraph # ',i);
    xlabel('set S_2 : red circles');
    ylabel('set S_1 : blue squares');
    title(S)
    set(gcf,'Position',[30+i*width,30,width,height]);
    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Scatter Plots


% wavelet coeffs of individual channels
w1 = f_w1(Colorednodes{1},1);
w2 = f_w1(Colorednodes{2},2);
w3 = f_w1(Colorednodes{3},3);
w4 = f_w1(Colorednodes{4},4);

% use same scale for high-pass coefficients
w_d = [w2(:); w4(:)];
w2 = abs(w2)/norm(w_d);
w3 = abs(w3)/norm(w_d);
w4 = abs(w4)/norm(w_d);

% axis([-97 -91 44 48])
% use same scale for original signal, low-pass wavelet coeffs and its low-pass reconsuction
% f_hat = sum(f_hat,2);
f_hat = f_hat(:,1);
Data = [f(:);f_hat(:,1)];
clim = [min(Data) max(Data)];
f = (f - clim(1))/(clim(2)-clim(1));
f = 2*f -1;
f_hat = (f_hat - clim(1))/(clim(2)-clim(1));
f_hat = 2*f_hat -1;
w_a = (w1 - clim(1))/(clim(2)-clim(1));
w_a = 2*w_a -1;
f_hat_d = abs(f_hat(:,2:end));
f_hat_d = f_hat_d /norm(f_hat_d(:));
% w_a = w1;
% set(gcf,'Position',[30,30+2*height+10,height,height]);
figure,
ti = get(gca,'TightInset');
set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
show_wavelet(f,xy(:,1),xy(:,2));
set(gcf,'Position',[30,(30+2*height),width,height]);
title('Input Signal')
axis auto

% wavelet coeffs
figure,
h = subplot(2,2,1);
set(h, 'pos',[0, 0.5 ,0.4,0.4]);
show_wavelet(w_a,xy(Colorednodes{1},1),xy(Colorednodes{1},2));
S = sprintf('%s%s%s', 'Wavelet coeffs of ',Channel_Name{1},' channel');
title(S)
axis auto

h = subplot(2,2,2);
set(h, 'pos',[0.5, 0.5 ,0.4,0.4]);
% set(gcf,'Position',[30+2*height,30+ height,height,height]);
show_wavelet(w2,xy(Colorednodes{2},1),xy(Colorednodes{2},2));
S = sprintf('%s%s', Channel_Name{2},' channel');
title(S)
axis auto

h = subplot(2,2,3);
set(h, 'pos',[0, 0 ,0.4,0.4]);
% set(gcf,'Position',[30+2*height,30+ height,height,height]);
show_wavelet(w3,xy(Colorednodes{3},1),xy(Colorednodes{3},2));
S = sprintf('%s%s', Channel_Name{2},' channel');
title(S)
axis auto

h = subplot(2,2,4);
set(h, 'pos',[0.5, 0 ,0.4,0.4]);
show_wavelet(w4,xy(Colorednodes{4},1),xy(Colorednodes{4},2));
S = sprintf('%s%s', Channel_Name{4},' channel');
title(S)
set(gcf,'Position',[30+width,(30+2*height),width,height]);
axis auto



%
% for i = 1:size(f_hat_d,2)
%     subplot(2,2,i+1)
%     S = sprintf('%s%s%s', 'Reconstructed signal from ',Channel_Name{i+1},' channel');
%     title(S)
%     set(gcf,'Position',[30+(i+2)*height,30,height,height]);
%     show_wavelet(f_hat_d(:,i),xy(:,1),xy(:,2));
%     axis auto
% end


% reconstrcuted signals
figure,
% h = subplot(2,2,1);
% set(h, 'pos',[0, 0.5 ,0.4, 0.4]);
show_wavelet(sum(f_hat,2),xy(:,1),xy(:,2));
S = sprintf('%s','Reconstruction from 1% detail coefficients');
title(S)
axis auto
%
% h = subplot(2,2,2);
% set(h, 'pos',[0.5, 0.5 ,0.4,0.4]);
% %     set(gcf,'Position',[30+(i+2)*height,30,height,height]);
% show_wavelet(f_hat_d(:,1),xy(:,1),xy(:,2));
% S = sprintf('%s%s', Channel_Name{2},' channel');
% title(S)
% axis auto
%
% h = subplot(2,2,3);
% set(h, 'pos',[0, 0 ,0.4,0.4]);
% %     set(gcf,'Position',[30+(i+2)*height,30,height,height]);
% show_wavelet(f_hat_d(:,2),xy(:,1),xy(:,2));
% S = sprintf('%s%s', 'HL',' channel');
% title(S)
% axis auto
%
% h = subplot(2,2,4);
% set(h, 'pos',[0.5, 0 ,0.4,0.4]);
% %     set(gcf,'Position',[30+(i+2)*height,30,height,height]);
% show_wavelet(f_hat_d(:,3),xy(:,1),xy(:,2));
% S = sprintf('%s%s',Channel_Name{4},' channel');
% title(S)
% set(gcf,'Position',[30+2.2*width,(30+2*height),2*width,2*height]);
% axis auto

% figure,
% set(gcf,'Position',[30+(i+1)*height,30,height,height]);
% show_wavelet(f_hat_d(:,2),xy(:,1),xy(:,2));
% axis auto

% for i = 1:Fmax
%     figure,
% %     set(gca,'ydir','reverse')
%    axis([-97 -91 44 48])
%     set(gcf,'Position',[30+i*height,30+2*height+10,height,height]);
%     mM = [min(f_hat(:,i)),max(f_hat(:,i))];
%     diffX = f_hat(:,i);
%     X1 = (diffX - mM(1))/(mM(2) - mM(1));
%     hold on
%     for indx = 1:length(C)
%         if all(C{indx}~=1)   % If at least one of the indices is 1,
%             % then it is an open region and we can't
%             % patch that.
%             patch(V(C{indx},1),V(C{indx},2),X1(indx)); % use color i.
%         end
%     end
%     axis([-97 -91 44 48])
%     colormap(gray)
%     S = sprintf('%s%s%s', 'Reconstructed signal from ',Channel_Name{i},' channel');
%     xlabel('Voronoi Tesselation Plot');
%     title(S)
% end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Voronoi Tesselation Plots
% [V,C] = voronoin(xy);
% figure,
% % voronoi(xy(:,1),xy(:,2));
% % set(gca,'ydir','reverse')
% axis([-97 -91 44 48])
% title('Input Signal')
% % set(gcf,'Position',[30,30+2*height+10,height,height]);
% set(gcf,'Position',[30+(i+1)*height,30,height,height]);
% mM = [min(f),max(f)];
% diffX = f;
% X1 = (diffX - mM(1))/(mM(2) - mM(1));
% hold on
% for indx = 1:length(C)
%     if all(C{indx}~=1)   % If at least one of the indices is 1,
%         % then it is an open region and we can't
%         % patch that.
%         patch(V(C{indx},1),V(C{indx},2),X1(indx)); % use color i.
%     end
% end
% axis([-97 -91 44 48])
% colormap(gray)
%
% for i = 1:Fmax
%     figure,
% %     set(gca,'ydir','reverse')
%    axis([-97 -91 44 48])
%     set(gcf,'Position',[30+i*height,30+2*height+10,height,height]);
%     mM = [min(f_hat(:,i)),max(f_hat(:,i))];
%     diffX = f_hat(:,i);
%     X1 = (diffX - mM(1))/(mM(2) - mM(1));
%     hold on
%     for indx = 1:length(C)
%         if all(C{indx}~=1)   % If at least one of the indices is 1,
%             % then it is an open region and we can't
%             % patch that.
%             patch(V(C{indx},1),V(C{indx},2),X1(indx)); % use color i.
%         end
%     end
%     axis([-97 -91 44 48])
%     colormap(gray)
%     S = sprintf('%s%s%s', 'Reconstructed signal from ',Channel_Name{i},' channel');
%     xlabel('Voronoi Tesselation Plot');
%     title(S)
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Delaunay triangulation plots
%
% minf_hat = min(f_hat(:));
% maxf_hat = max(f_hat(:));
% minf_w = min(f_w(:));
% maxf_w = max(f_w(:));
% minf = min([minf_hat minf_w]);
% maxf = max([maxf_hat maxf_w]);
% figure,
% gridDelaunay = delaunay(xy(:,1),xy(:,2));
% trimesh(gridDelaunay,xy(:,1),xy(:,2),f)
% view(0,90)
% colormap('cool')
% v = axis;
% caxis([minf maxf])
% % axis([-98 -88 43 50])
% colorbar
% set(gcf,'Position',[30,30+3*height+20,height,height]);
% title(' Graph-signal')
% xlabel('Delaunay triangulation plot')
% axis(v)
% % f_w = abs(f_w);
% % maxf = max(f_w(:));
% % minf = min(f_w(:));
% % f_w = (f_w - minf)/(maxf - minf);
% warning off
% % for i = 1:Fmax
% %     if length(Colorednodes{i}) >2 % need atleast 3 points for triangulation
% %         figure,
% %         gridDelaunay = delaunay(xy(Colorednodes{i},1),xy(Colorednodes{i},2));
% %         trimesh(gridDelaunay,xy(Colorednodes{i},1),xy(Colorednodes{i},2),f_w(Colorednodes{i},i))
% %         view(0,90)
% %         colormap('cool')
% %         set(gcf,'Position',[30+i*height,30+3*height+20,height,height]);
% %         S = sprintf('%s%s%s',Channel_Name{i},' channel ', 'wavelet coefficients');
% %         xlabel('Delaunay triangulation plot');
% %         title(S)
% %         caxis([minf maxf])
% %         colorbar
% %         axis(v);
% %     end
% % end
%
% gridDelaunay = delaunay(xy(:,1),xy(:,2));
%
% for i = 1:Fmax
%     figure,
%     trimesh(gridDelaunay,xy(:,1),xy(:,2),f_hat(:,i))
%     view(0,90)
%     colormap('cool')
%     set(gcf,'Position',[30+i*height,30+3*height+20,height,height]);
%     S = sprintf('%s%s%s', 'Reconstructed signal from ',Channel_Name{i},' channel');
%     xlabel('Delaunay triangulation plot');
%     title(S)
%     caxis([minf maxf])
%     colorbar
%     axis(v);
% end

