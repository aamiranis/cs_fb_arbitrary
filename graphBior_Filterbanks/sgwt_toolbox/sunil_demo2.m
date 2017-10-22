% demo 4
n = 16;
m = n;
N = n^2;
% find dct basis functions
% B1 = zerosh
% A1 = square_grid(n);
A1 = sgwt_meshmat(n);
% D1 = diag(sum(A1,2));
L1 = 4*eye(256) - A1;
[U1 Lam1] = eig(L1);
Lam1 = diag(Lam1);
[Lam1 I] = sort(Lam1);
U1 = U1(:,I);
plot(U1(:,3))
% 
% f = rand(256,1);%1:256; f = f';
% I = reshape(f,16,16);
% % I = (1:16)';
% % I = repmat(I,1,16);
% % f = I(:);
% alpha1 = f;
% alpha2 = dct2(I);
% alpha2 = alpha2(:);
% figure,
% plot(sort(alpha1));
% hold on
% plot(sort(alpha2),'r');
% 


% A = sgwt_meshmat(n);
% D = ones(N,1);
% for i = 1:m
%     A(i,m+i) = 2;
%     D(i) = sqrt(2);
%     A((n-1)*m+i, (n-2)*m+i) = 2;
%     D((n-1)*m+i) = sqrt(2);
% end
% for j = 1:n
%     A((j-1)*m+1,(j-1)*m+2) = 2;
%     D((j-1)*m+1) = sqrt(2);
%     A(j*m,j*m-1) = 2;
%     D(j*m) = sqrt(2);
% end
% L2 = 4*eye(N) - A;
% [U2 Lam2] = eig(L2);
% Lam2 = real(diag(Lam2));
% [Lam2 I] = sort(Lam2);
% U2 = real(U2(:,I));
% 
% D(1) = 2;
% D(m) = 2;
% D((n-1)*m+1) = 2;
% D(N) = 2;
% D = spdiags(D,0,N,N);
% A = D^(-1)*A*D;
% 
% L3 = 4*eye(N) - A;
% [U3 Lam3] = eig(L3);
% Lam3 = diag(Lam3);
% [Lam3 I] = sort(Lam3);
% U3 = U3(:,I);
% 
% figure(1),
% plot(Lam1,'b*-')
% hold on
% % plot(Lam2,'rs-')
% plot(Lam3,'mp-')
% legend('original spectrum','symmetrized spectrum')
% figure(2)
% imagesc(U1'*U3);
% title('original vs extended');
% 
% % figure(3)
% % imagesc(U2'*U3)
% % title('extended vs symmetrized')
