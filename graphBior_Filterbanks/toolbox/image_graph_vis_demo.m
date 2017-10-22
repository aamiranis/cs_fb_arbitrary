% Image graph visualization demo
clear 
clc
s_im = 8;
N = s_im*s_im;
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
Ar=sparse([cind1,nind1],[nind1,cind1],ones(1,2*numel(ni)),N,N);
Ar = full(Ar);
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
Ad=sparse([cind2,nind2],[nind2,cind2],ones(1,2*numel(ni)),N,N);
Ad = full(Ad);

% vertical
ci=[alli];%alli];
cj=[allj];%allj];
ni=[alli];%  ; alli+1];
nj=[allj+1];%; allj];
% prune edges at boundary
valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2));
ni=ni(valid);
nj=nj(valid);
ci=ci(valid);
cj=cj(valid);
cind3=dim(1)*(cj-1)+ci;
nind3=dim(1)*(nj-1)+ni;
Av=sparse([cind3,nind3],[nind3,cind3],ones(1,2*numel(ni)),N,N);
Av = full(Av);

% horizontal 
Ah = Ar - Av;
xy = [alli,allj]/(2*N);




% 8-connected coloring
beta(:,1) = image_downsampling_fn(s_im,'rectangle');
beta(:,2) = image_downsampling_fn(s_im,'diamond');
beta_temp = 0.5*(1+beta);
F = 2*beta_temp(:,1 ) + beta_temp(:,2);
clear beta_temp
theta = 2;
F = 2^theta - F ;
Fmax = 2^theta;
figure,
axis tight
gplot((Ar+Ad),xy)
hold on
S1 = find(F == 1);
S2 = find(F == 2);
S3 = find(F == 3);
S4 = find(F == 4);
plot(xy(S1,1),xy(S1,2),'ks','MarkerFaceColor','b','MarkerSize',10);
plot(xy(S2,1),xy(S2,2),'ks', 'MarkerFaceColor','k','MarkerSize',10);
plot(xy(S3,1),xy(S3,2),'ks', 'MarkerFaceColor','r','MarkerSize',10);
plot(xy(S4,1),xy(S4,2),'ks', 'MarkerFaceColor','g','MarkerSize',10);

figure,
axis tight
gplot(Ar,xy)
hold on
S1 = find(F == 1);
S2 = find(F == 2);
S3 = find(F == 3);
S4 = find(F == 4);
plot(xy(S1,1),xy(S1,2),'ks','MarkerFaceColor','b','MarkerSize',10);
plot(xy(S2,1),xy(S2,2),'ks', 'MarkerFaceColor','b','MarkerSize',10);
plot(xy(S3,1),xy(S3,2),'ks', 'MarkerFaceColor','r','MarkerSize',10);
plot(xy(S4,1),xy(S4,2),'ks', 'MarkerFaceColor','r','MarkerSize',10);


figure,
axis tight
gplot(Ad,xy)
hold on
S1 = find(F == 1);
S2 = find(F == 2);
S3 = find(F == 3);
S4 = find(F == 4);
plot(xy(S1,1),xy(S1,2),'ks','MarkerFaceColor','b','MarkerSize',10);
plot(xy(S2,1),xy(S2,2),'ks', 'MarkerFaceColor','r','MarkerSize',10);
plot(xy(S3,1),xy(S3,2),'ks', 'MarkerFaceColor','b','MarkerSize',10);
plot(xy(S4,1),xy(S4,2),'ks', 'MarkerFaceColor','r','MarkerSize',10);


figure,
axis tight
gplot(Av,xy)
hold on
S1 = find(F == 1);
S2 = find(F == 2);
S3 = find(F == 3);
S4 = find(F == 4);
plot(xy(S1,1),xy(S1,2),'ks','MarkerFaceColor','b','MarkerSize',10);
plot(xy(S2,1),xy(S2,2),'ks', 'MarkerFaceColor','r','MarkerSize',10);
plot(xy(S3,1),xy(S3,2),'ks', 'MarkerFaceColor','b','MarkerSize',10);
plot(xy(S4,1),xy(S4,2),'ks', 'MarkerFaceColor','r','MarkerSize',10);

figure,
axis tight
gplot(Av,xy)
hold on
S1 = find(F == 1);
S2 = find(F == 2);
S3 = find(F == 3);
S4 = find(F == 4);
plot(xy(S1,1),xy(S1,2),'ks','MarkerFaceColor','b');
plot(xy(S2,1),xy(S2,2),'ks', 'MarkerFaceColor','r');
plot(xy(S3,1),xy(S3,2),'ks', 'MarkerFaceColor','b');
plot(xy(S4,1),xy(S4,2),'ks', 'MarkerFaceColor','r');


