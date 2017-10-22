function output = non_linear_approx_standard(varargin)
clc
% close all
disp('This function computes wavelet coefficients of any Image using a 2-dim separable standard CDF 9/7 wavelet transform.')
addpath Graph_Generators/
addpath Graph_kernels/
addpath sgwt_toolbox/
addpath toolbox/
addpath exportfig/

p = inputParser;
p.addOptional('filename','ballet.png',@(f) exist(f,'file'));
p.addParamValue('maxLevel',4,@isscalar);
p.parse(varargin{:});
opt = p.Results;


% Set parameters
filename = opt.filename;
tempI = imfinfo(filename);
filetype = tempI.Format;
max_level = opt.maxLevel; % number of decomposition levels

%% Graph Signal
Data = imread(filename,filetype);
if length(size(Data)) == 3
    Data = rgb2gray(Data);
end
Data = double(Data);%/255;
[m n] = size(Data);
m = min([m,n]);
s_im = floor(m/2^max_level)*2^max_level;
% s_im = 512;
% s_im = 512;
% s_im = 10;
Data = Data(1:s_im,1:s_im);
% Data = Data/norm(Data(:));
Data_w = wavecdf97(Data,max_level);
N(1) = s_im*s_im;
for level = 2:max_level
    N(level) = N(level-1)/4;
end
dim = sqrt(N(max_level))/2;
output{1} = Data_w;
tempw = Data_w(1:dim,(dim+1):s_im);
tempw1 = Data_w((dim+1):s_im,:);
tempw = [tempw(:); tempw1(:)];
% hist(abs(tempw),64);
% set(gcf, 'Color', 'w');
% V = axis;
% V(4) = 2*10^5;
% V(2) = 9*10^-3;
% axis(V);
% outpfile = sprintf('%s%s%s','Results/',filename,'_hist_standardCDF.png');
% export_fig(outpfile,'-png','-m2');
norm0 = length(find(abs(tempw)>10^-6));%/length(tempw);
output{2} = [norm0 norm(tempw,1) norm(tempw,2) norm(tempw,Inf) norm(tempw,-Inf)];
output{3} = abs(tempw);
% tempw = Data_w;
% tempw(1:dim,1:dim) = 0;
% Data_w = Data_w - tempw;
% Data_w = rescale(Data_w) + rescale(abs(tempw));
% % Data_w = Data_w.^2;
% outpfile = sprintf('%s%s%s%s%d%s','Results/',filename,'_standardCDF_wavelet_coeffs','_levels',max_level,'.png');
% imwrite(Data_w,outpfile,'png');
% temp_w = Data_w1;
% len = s_im/2^max_level;
% temp_w(1:len,1:len) = 0;
% Data_w = Data_w1 - temp_w;
% Data_w = rescale(Data_w) + rescale(abs(100*temp_w));
% figure,imagesc(Data_w);
% colormap(gray)
origImg = Data;
f = Data(:);
nnz = 0:0.01:0.18;
MSE = zeros(1,length(nnz));
SSIM = zeros(1,length(nnz));
PSNR = zeros(1,length(nnz));
SNR = zeros(1,length(nnz));
Data_w1 = Data_w;
for iter = 1:length(nnz)
    nnz_factor = nnz(iter);
    Data_w = Data_w1;
    % save ballet_image_CDF97 Data_w
    temp_w = Data_w;
    len = s_im/2^max_level;
    temp_w(1:len,1:len) = 0;
    coeffs = temp_w(:);
    thresh = sort(abs(coeffs),'descend');
    %     nCoeffs = nnz_factor -1;
    nCoeffs = s_im^2 - len^2;
    nCoeffs = floor(nCoeffs*nnz_factor);
    if nnz_factor ==0
        coeffs = 0*coeffs;
    else
        thresh = thresh(nCoeffs+1);
        coeffs(abs(coeffs) <thresh) = 0;
    end
    temp_w = reshape(coeffs,s_im,s_im);
    temp_w(1:len,1:len) = Data_w(1:len,1:len);
    Data_w = temp_w;
    Data_hat = wavecdf97(Data_w,-max_level);
    if iter ==4
         low_in = min(Data_hat(:));
        high_in = max(Data_hat(:));
        high_out = 255;
        Data_hat1 = high_out*(Data_hat - low_in)/(high_in - low_in);
        figure,imagesc(Data_hat1)
%         axis([95 160 85 160]);
        colormap(gray)
        colorbar
        set(gca,'xTick',[])
        set(gca,'yTick',[])
        title(sprintf('%s%1.2f%s','Reconstruction using ',nnz_factor ,' highpass coefficients using standard CDF'));
    end
    distImg = Data_hat;
    f_hat1 = distImg(:);
    MSE(iter) = MeanSquareError(origImg, distImg);
    % disp('Mean Square Error = ');
    % disp(MSE);
    
    %Peak Signal to Noise Ratio
    PSNR(iter) = PeakSignaltoNoiseRatio(origImg, distImg);
    % disp('Peak Signal to Noise Ratio = ');
    % disp(PSNR);
    
    SSIM(iter) = ssim(origImg,distImg);
    MSE1 = (f - f_hat1).^2;
    MSE1 = sum(MSE1);
    SNR(iter) = 10*log10(norm(f)^2./MSE1);
end
output{4} = struct('SNR',SNR,'PSNR',PSNR, 'SSIM',SSIM,'MSE', MSE);