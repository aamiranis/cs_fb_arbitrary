LN = 20;

% A = square_grid(N);
% Q=load('minnesota.mat');
% xy=Q.xy;
% A=Q.A;
% A = full(A);
% random graph
clear all
close all
clc
% N = 400;
% p = 0.05;
% A = rand(N);
% G = A < p;
% G = triu(G,1);
% G = G + G';
% Sum = sum(G,2);
% ind = find(Sum ~=0);
% A = G(ind,ind); 
% A = line_graph(400);
A = square_grid(20);
N = length(A);
D = diag(sum(A,2));
L = D - A;
[U Lam] = eig(L);
Lam = diag(Lam);
% 
% Lam1 = Lam;
% Lam1(2:N) = 1./Lam(2:N);
% Lplus = U*diag(Lam1)*U';
% Dist = 0*L;
% for i = 1: N
%     Dist(i,:) = Lplus(i,i) + diag(Lplus)' - 2*Lplus(i,:);
% end
% Lnorm = D^-(0.5)*L*D^(-0.5);
% [U Lam] = eig(Lnorm);
Lam1 = Lam;
Lam1(2:N) = 1./Lam(2:N);
Lplus1 = U*diag(Lam1)*U';
Lplus = pinv(L);
Lplus = 0.5*(Lplus  + Lplus');
Dist = 0*L;
for i = 1: N
    Dist(i,:) = Lplus(i,i) + diag(Lplus)' - 2*Lplus(i,:);
end
Dist = 0.5*(Dist+Dist');
T = 0*L;
dmax = max(Dist(:));
for i= 1:N
    T(i,:) = max_hat(Dist(i,:),dmax);
    T(i,i) = T(i,i)-sum(T(i,:));
%     pos = find(T(i,:)>0);
%     neg= find(T(i,:)<0);
%     sum_pos = sum(T(i,pos));
%     sum_neg = -sum(T(i,neg));
%     T(i,pos) = T(i,pos)/sum_pos;
%     T(i,neg) = T(i,neg)/sum_neg; 
end
% T = 0.5*(T+T');
% Dist1 = zeros(400,40);
% for i = 1: 40
%     P = inv(D)*A;
%     Pi = P^i;
%     DPi = Pi*inv(D);
%     Dist1(:,i)= DPi(1,:)';
% end

    
    
% for i = 1:N
%     T(i,Dist(i,:) < 0.5) = 0.5;
%     T(i,(Dist(i,:) >= 0.5 & Dist(i,:) < 1) ) = -0.5;
% end
% for i = 1:N
% T(i,:) = exp(-0.5*Dist(i,:)).* cos(Dist(i,:));
% end
% Dist1 = pinv(Dist);
% for i = 1:N
% T(i,:) = Dist1(i,:);
% end
% T = Dist + Dist.*Dist + Dist.*Dist.*Dist;
G = U'*T*U;
for i = 1:N
    plot(abs(G(i,:)))
    pause
end
% for i = 1: 400
%     plot(Dist1(i,:),'*')
%     hold on
% end

% Vg = sum(A(:));
% Dist = Vg*Dist;
% A1 = A;
% 
% 
% Atemp = A + A^2;
% A2 = double(Atemp >0);
% A2 = A2-A1; % A2 is a ring adjacency matrix
% 
% Atemp = Atemp + A^3;
% A3 = double(Atemp>0);
% A3 = A3 - A2 - A1;
% 
% for i = 1:N
%     n1 = find(A1(i,:));
%     n2 = find(A2(i,:));
%     n3 = find(A3(i,:));
%     n4 = setdiff(1:N,[n1,n2,n3]);
%     plot(Dist(i,[n1,n2,n3,n4]),'*')
%     pause
% end
    
    

% [U1 Lam1]= eig(Dist);
% [x y] = draw_dot(A);
% for i = 1:N
%     dist = Dist(i,:)';
%     dist = reshape(dist,20,20);
%     imagesc(dist);
%     pause
% end

