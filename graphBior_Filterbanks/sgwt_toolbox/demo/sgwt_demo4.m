% sgwt_demo3 : Image decomposition with SGWT wavelets based on local adjacency.
%
% This demo builds the SGWT transform on a graph representing
% adjacency on a pixel mesh with 4-nearest neighbor connectivity.
% This demonstrates inverse on problem with large dimension.
%
% The demo loads an image file and decomposes the image with the SGWT,
% showing the coefficients as images at each scale. The demo does not show
% the individual wavelets (this could be done by replacing the input
% image by a "delta image" with a single unit nonzero pixel) .
%
% The inverse is then computed, from the original coefficients as well as
% from a modified set of coefficients where only coefficients at one
% scale are preserved. This shows that the SGWT can generate a
% multiresolution decomposition for images. We don't claim that this
% particular local-adjacency based transform is better for image
% processing than other available wavelet image decompositions, but it
% demonstrates the flexibility of the SGWT.

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond.
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

% function sgwt_demo4
close all;
fprintf('Welcome to SGWT demo #4\n');
% % load image
% imname='squares.jpg';
% fprintf('loading image %s\n',imname);
% im = double( imread(imname));
% [s1 s2] = size(im);
% s_im = min(s1,s2);
s_im = 300;
im = zeros(s_im);
for i = 1:1:s_im
    im(:,i) = cos(pi*i/4);
end
% for i = 1:1:s_im
%     im(i,:) = im(i,:)*cos(pi*i/4);
% %     im(i,:) = cos(pi*i/s_im);
% end
% im1 = imrotate(im,-45,'bilinear','crop');
% diff = ceil(0.5*(s_im *sqrt(2) - s_im));
% x = diff:(s_im-diff);
% % s_im = 100;%floor(s_im/2)*2;
% % im = im(1:s_im,1:s_im)';
% im = im1(x,x);
s_im = length(im);     

% im = zeros(s_im);   
% im(s_im/2,:) = 255;
% im = 255*(im>0);
% build mesh adjacency graph
fprintf('Building mesh adjacency graph\n');
% boundary = 'rectangle';
% boundary = 'torus';
% boundary = 'diamond';
boundary = 'horizontal';
A=sgwt_meshmat(size(im),'boundary',boundary);

%k-regularize matrix A
% A = kregularize(A);
A = 0.5*(A+A');
% transform
fprintf('Calculating graph Laplacian\n');
L=sgwt_laplacian(A);

% L1=sgwt_laplacian(A);
% L = 4*speye(size(A)) - A;
fprintf('Measuring largest eigenvalue, lmax = ');
lmax=sgwt_rough_lmax(L);
% lmax = 4;
arange=[0,lmax];
fprintf('%g\n',lmax);

Nscales=2;
fprintf('Designing transform in spectral domain\n');
% [g,gp,t]=sgwt_filter_design(lmax,Nscales);
g{1} = @(x)heaviside(lmax/2 - x); % low-pass kernel
g{2} = @(x)heaviside(x - lmax/2); % high-pass kernel
m=4*s_im; % order of polynomial approximation
fprintf('Computing Chebyshev polynomials of order %g for fast transform \n',m);
for k=1:numel(g)
    c{k}=sgwt_cheby_coeff(g{k},m,m+1,arange);
end

fprintf('Computing forward transform\n');
wpall=sgwt_cheby_op(im(:),L,c,arange);

% L = full(L);
% [U Lam] = eig(L);
% g_low = [ones(1,s_im^2/2),  zeros(1,s_im^2/2)];
% % g_low = rand(1,s_im^2).*g_low;
% % g_high = fliplr(g_low);
% % T_low = U*diag(g_low)*U';
% f = U*g_low';
% immm = reshape(f,size(im));
% T_high = U*diag(g_high)*U';
% wpall{1} = T_low*im(:);
% wpall{2} = T_high*im(:);
E = [];
O= [];
if strcmp(boundary,'rectangle')
    %     E = randsample(s_im^2,s_im^2/2);
    %     O = setdiff(1:s_im^2,E);
    for j = 1:2:s_im
        E = [E ((j-1)*s_im + (1:2:s_im))];
        O = [O ((j-1)*s_im + (2:2:s_im))];
    end
    
    for j = 2:2:s_im
        E = [E ((j-1)*s_im + (2:2:s_im))];
        O = [O ((j-1)*s_im + (1:2:s_im))];
    end
elseif strcmp(boundary,'diamond')
    for j = 1:2:s_im
        E = [E ((j-1)*s_im + (1:1:s_im))];
        O = [O (j*s_im + (1:1:s_im))];
    end
elseif strcmp(boundary,'horizontal')
    for j = 1:2:s_im
        E = [E ((j-1)*s_im + (1:1:s_im))];
        O = [O (j*s_im + (1:1:s_im))];
    end
end
wpall{1}(O) = 0;
wpall{2}(O) = 0;



fprintf('Computing inverse transform\n');

im_hat_low=sgwt_cheby_op(wpall{1},L,c{1},arange);
im_hat_high= sgwt_cheby_op(wpall{2},L,c{2},arange);
% imm_low = reshape(wpall{1},size(im));
% imm_high = reshape(wpall{2},size(im));
imr_low = reshape(im_hat_low,size(im));
imr_high = reshape(im_hat_high,size(im));
% MSE = im - imr_low;
% MSE = MSE.^2;
% MSE = sum(MSE(:))/s_im^2;
% % invert with all subbands
% fprintf('Computing inverse transform with all coefficients\n');
% imr1=sgwt_inverse(wpall,L,c,arange);
% imr1=reshape(imr1,size(im));
%
% ks=2; % scale at which to keep coefficients, set all others to zero.
% fprintf('\nsetting all coefficients to zero except wavelet scale %g\n',ks-1);
% % invert with only one scale
% for k=1:numel(wpall)
%     wpall2{k}=zeros(size(wpall{k}));
% end
% wpall2{ks}=wpall{ks};
% fprintf('Computing inverse transform with coefficients from wavelet scale %g only\n',ks-1);
% imr2=sgwt_inverse(wpall2,L,c,arange);
% imr2=reshape(imr2,size(im));
%
% %% display results
% figure(1)
% set(gcf,'position',[ 5   730   350   350]);
% sgwt_show_im(im)
% title('original image');
% set(gcf,'menubar','none')
% figure(2)
% set(gcf,'position',[365 730 350 350]);
% sgwt_show_im(imr1)
% title('reconstuction from all coefficients');
% set(gcf,'menubar','none')
%
% figure(3)
% set(gcf,'position',[725 730 350 350]);
% sgwt_show_im(imr2);
% title(sprintf('reconstruction only from wavelets at scale %g',ks-1));
% set(gcf,'menubar','none')
%
% figure(4)
% set(gcf,'position',[0 0 1150 700]);
% set(gcf,'menubar','none')
% for k=1:Nscales
%     subplot(1,2,k);
%     sgwt_show_im(reshape(wpall{k},size(im)));
%     if k==1
%         title('Scaling function coefficients');
%     else
%         title(sprintf('Wavelet coefficients at scale %g',k-1));
%     end
% end
% ks = 2;
% %% display results
% figure(1)
% set(gcf,'position',[ 5   730   350   350]);
% sgwt_show_im(im)
% title('original image');
% % set(gcf,'menubar','none')
% 
% figure(2)
% set(gcf,'position',[725 730 350 350]);
% sgwt_show_im(imr_high);
% title(sprintf('reconstruction only from wavelets at scale %g',ks-1));
% % set(gcf,'menubar','none')
% 
% figure(3)
% set(gcf,'position',[365 730 350 350]);
% % set(gcf,'menubar','none')
% sgwt_show_im(imr_low);
% title('Scaling function coefficients');
% 
% lambda=linspace(0,lmax,1e3);
% N = length(lambda);
% h1 = [ones(1,N/2) zeros(1,N/2)];
% h2 = [zeros(1,N/2) ones(1,N/2)];
% h1_hat = sgwt_cheby_eval(lambda,c{1},arange);
% h2_hat = sgwt_cheby_eval(lambda,c{2},arange);
% figure(5)
% % set(gcf,'position',[425,580,600,250])
% plot(lambda,h1,lambda,h1_hat);
% legend('Exact Scaling kernel','Chebyshev polynomial approximation');
% 
% figure(6)
% % set(gcf,'position',[425,580,1150,700])
% plot(lambda,h2,lambda,h2_hat);
% legend('Exact Highpass kernel','Chebyshev polynomial approximation');
% 
% 
