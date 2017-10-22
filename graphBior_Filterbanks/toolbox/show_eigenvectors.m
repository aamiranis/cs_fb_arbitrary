% eigenvector demo
clear all
clc
close all
[A xy f] = uniform_graphs();  % uniformly distributed graphs
scrsz = get(0,'ScreenSize');
height = scrsz(4)/4;
figure,
gplot(A,xy,'r.-')
title('Original Graph')
set(gcf,'Position',[30,30,height,height]);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

% compute eigenvectors
D = diag(sum(A,2));
L = D - A;
L = 0.5*(L + L');
[U Lam] = eig(L);
Lam = diag(Lam);
ind = [10 40 300 500];
figure,
for i = 1:length(ind)
    subplot(2,2,i)
    gplot(A,xy,'r.-')
    hold on
    u = U(:,ind(i));
    u = double(u > 0);
    u = 2*u -1;
     show_wavelet(u,xy(:,1),xy(:,2))
    colormap(gray)
end

    
   
