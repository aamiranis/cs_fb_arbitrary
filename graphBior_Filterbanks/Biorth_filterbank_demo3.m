% L-1 norm comparisons of all graphBior designs with standard wavelets
close all
% filename = 'Lichtenstein.png';
filename = 'coins.png';
tempI = imfinfo(filename);
filetype = tempI.Format;
% graphBior asym GC off edgemap 0
% output =  non_linear_approx_graphBior(filename,'LType','asym','GC','off','edgemap',false);
% wav_coeffs1 = output{1};
% wav_norm1 = output{2};
% dist1 = output{3};
% SNRdata1 = output{4};
% % graphBior asym GC off edgemap 1
% output = non_linear_approx_graphBior(filename,'LType','asym','GC','off','edgemap',true);
% wav_coeffs2 = output{1};
% wav_norm2 = output{2};
% dist2 = output{3};
% SNRdata2 = output{4};

% graphBior asym GC on edgemap 0
output =  non_linear_approx_graphBior(filename,'LType','asym','GC','on','edgemap',false,'maxLevel',4);
wav_coeffs3 = output{1};
wav_norm3 = output{2};
dist3 = output{3};
SNRdata3 = output{4};

% graphBior asym GC on edgemap 1
output =  non_linear_approx_graphBior(filename,'LType','asym','GC','on','edgemap',true,'maxLevel',4);
wav_coeffs4 = output{1};
wav_norm4 = output{2};
dist4 = output{3};
SNRdata4 = output{4};

output = non_linear_approx_standard(filename,'maxLevel',4);
wav_coeffs5 = output{1};
wav_norm5 = output{2};
dist5 = output{3};
SNRdata5 = output{4};

output =  non_linear_approx_graphBior(filename,'LType','sym','GC','on','edgemap',false,'maxLevel',4);
wav_coeffs6 = output{1};
wav_norm6 = output{2};
dist6 = output{3};
SNRdata6 = output{4};

output =  non_linear_approx_graphBior(filename,'LType','sym','GC','on','edgemap',true,'maxLevel',4);
wav_coeffs7 = output{1};
wav_norm7 = output{2};
dist7 = output{3};
SNRdata7 = output{4};

max_level = 4;
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
low_in = min(f);
high_in = max(f);
high_out = 255;
distImg1 = high_out*(double(Data) - low_in)/(high_in - low_in);
figure,
imagesc(distImg1);
% axis([95 160 85 160]);
colorbar
colormap(gray)
set(gca,'xTick',[])
set(gca,'yTick',[])
title('Original Image')
figure,
nnz = 0:0.01:0.18;
% plot(nnz,SNRdata1.PSNR,'.-')
hold on
% plot(nnz,SNRdata2.PSNR,'r.-')
plot(nnz,SNRdata3.PSNR,'k.-')
plot(nnz,SNRdata4.PSNR,'m.-')
plot(nnz,SNRdata5.PSNR,'c.-')
plot(nnz,SNRdata6.PSNR,'r^-')
plot(nnz,SNRdata7.PSNR,'^-')
% legend('asym\_GCoff,edge\_off','asym\_GCoff,edge\_on','asym\_GCon,edge\_off','asym\_GCon,edge\_on','CDF9/7','sym\_GCoff,edge\_on','sym\_GCon,edge\_on');
legend('zeroDC graphBior','edge-aware zeroDC graphBior','CDF9/7','nonzeroDC graphBior','edge-aware nonzeroDC graphBior');
% legend('asym\_GCoff,edge\_off','asym\_GCon,edge\_off','CDF9/7','sym\_GCoff,edge\_on','sym\_GCon,edge\_on');
title('PSNR Comparison')
xlabel('fraction of detail coeffs used');
ylabel('PSNR(dB)');

figure,
nnz = 0:0.01:0.18;
% plot(nnz,SNRdata1.SSIM,'.-')
hold on
% plot(nnz,SNRdata2.SSIM,'r.-')
plot(nnz,SNRdata3.SSIM,'k.-')
plot(nnz,SNRdata4.SSIM,'m.-')
plot(nnz,SNRdata5.SSIM,'c.-')
plot(nnz,SNRdata6.SSIM,'r^-')
plot(nnz,SNRdata7.SSIM,'^-')
% legend('asym\_GCoff,edge\_off','asym\_GCoff,edge\_on','asym\_GCon,edge\_off','asym\_GCon,edge\_on','CDF9/7','sym\_GCoff,edge\_on','sym\_GCon,edge\_on');
legend('zeroDC graphBior','edge-aware zeroDC graphBior','CDF9/7','nonzeroDC graphBior','edge-aware nonzeroDC graphBior');
% legend('asym\_GCoff,edge\_off','asym\_GCon,edge\_off','CDF9/7','sym\_GCoff,edge\_on','sym\_GCon,edge\_on');
title('SSIM Comparison')
xlabel('fraction of detail coeffs used');
ylabel('SSIM');
