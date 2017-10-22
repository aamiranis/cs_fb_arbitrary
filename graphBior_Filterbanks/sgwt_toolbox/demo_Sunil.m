% function demo_Sunil
close all
boundary = 'rectangle';
Ar=sgwt_meshmat(10,'boundary',boundary);
boundary = 'diamond';
Ad=sgwt_meshmat(10,'boundary',boundary);
s_im = sqrt(length(Ar));
E1 = [];
O1= [];
E2 = [];
O2= [];

for j = 1:2:s_im
    E1 = [E1 ((j-1)*s_im + (1:2:s_im))];
    O1 = [O1 ((j-1)*s_im + (2:2:s_im))];
end

for j = 2:2:s_im
    E1 = [E1 ((j-1)*s_im + (2:2:s_im))];
    O1 = [O1 ((j-1)*s_im + (1:2:s_im))];
end

for j = 1:2:s_im
    E2 = [E2 ((j-1)*s_im + (1:1:s_im))];
    O2 = [O2 (j*s_im + (1:1:s_im))];
end

[x y] = draw_dot(Ar);
gplot(Ar,[x' y'])
hold on
% plot(x(E1),y(E1),'bs','MarkerFaceColor','b','MarkerSize',10);
% plot(x(O1),y(O1),'ro','MarkerFaceColor','r','MarkerSize',10);
title('rectangular graph')

Arh = Ar;
for i = 1:s_im^2-1
    Arh(i,i+1) = 0;
    Arh(i+1,i) = 0;
end
figure,gplot(Arh,[x' y'])
hold on
plot(x(E2),y(E2),'bs','MarkerFaceColor','b','MarkerSize',10);
% plot(x(O2),y(O2),'b.','MarkerFaceColor','r','MarkerSize',10);
title('Horizontal Decomposition')

Arv = Ar - Arh;
figure,gplot(Arv,[x' y'])
hold on
plot(x(E2),y(E2),'bs','MarkerFaceColor','b','MarkerSize',10);
% plot(x(O2),y(O2),'b.','MarkerFaceColor','r','MarkerSize',10);
title('Vertical Decomposition')




% [x y] = draw_dot(Ad);
figure,gplot(Ad,[x' y'])
hold on
plot(x(E2),y(E2),'bs','MarkerFaceColor','b','MarkerSize',10);
plot(x(O2),y(O2),'ro','MarkerFaceColor','r','MarkerSize',10);
title('diamond graph')
dim = [s_im,s_im];
[alli,allj]=find(ones(s_im));
ci=alli;
cj=allj;
ni=alli+1;
nj=allj+1;
% prune edges at boundary
valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2));
ni=ni(valid);
nj=nj(valid);
ci=ci(valid);
cj=cj(valid);
cind=dim(1)*(cj-1)+ci;
nind=dim(1)*(nj-1)+ni;
Ad1=sparse([cind,nind],[nind,cind],ones(1,2*numel(ni)),s_im^2,s_im^2);
figure,gplot(Ad1,[x' y'])
hold on
plot(x(E1),y(E1),'bs','MarkerFaceColor','b','MarkerSize',10);
plot(x(O1),y(O1),'ro','MarkerFaceColor','r','MarkerSize',10);
title('main diag decomposition')
Ad2 = Ad - Ad1;
figure,gplot(Ad2,[x' y'])
hold on
plot(x(E1),y(E1),'bs','MarkerFaceColor','b','MarkerSize',10);
plot(x(O1),y(O1),'ro','MarkerFaceColor','r','MarkerSize',10);
title('off diag decomposition')

A = Ad+ Ar;
[x y] =draw_dot(A);
gplot(A,[x' y'])
title('8-connected graph')

%
% % Downsampling
% Ar1 = Ar(E1,O1);
% Ar1 = Ar1*Ar1';
% Ar1 = Ar1 - diag(diag(Ar1));
% [x y] = draw_dot(Ar1);
% gplot(Ar1,[x' y'])
% title('rectangular graph after downsampling')
%
% Ad1 = Ad(E2,O2);
% Ad1 = Ad1*Ad1';
% Ad1 = Ad1 - diag(diag(Ad1));
% [x y] = draw_dot(Ad1);
% gplot(Ad1,[x' y'])
% title('diamond graph after downsampling')
%




% Ad = kextension(Ad);
% Ar = kextension(Ar);
% draw_dot(Ar);
% figure,draw_dot(Ad);



%     function B = kextension(A)
%         n = length(A);
%         D = full(sum(A,2));
%         extra = 4 - D;
%         extra = sum(extra); % require these many extra nodes
%         B = zeros(n+extra);
%         B(1:n,1:n) = A;
%         M = cumsum(4-D);
%         m0 = 1;
%         for i = 1:n
%             B(i,n + (m0:M(i))) = 1;
%             B(n + (m0:M(i)),i) = 1;
%             m0 = M(i)+1;
%         end
%     end


% end