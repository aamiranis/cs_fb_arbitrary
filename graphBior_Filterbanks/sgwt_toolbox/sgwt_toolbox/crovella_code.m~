%% filter-designs based on various methods

% create graph 
N = 20;

% A = square_grid(N);
% Q=load('minnesota.mat');
% xy=Q.xy;
% A=Q.A;
% A = full(A);
% random graph
clear all
close all
clc
N = 400;
p = 0.05;
A = rand(N);
G = A < p;
G = triu(G,1);
G = G + G';
Sum = sum(G,2);
ind = find(Sum ~=0);
A = G(ind,ind); 
N = length(A);
D = diag(sum(A,2));
L = D - A;
[U Lam] = eig(L);
Lam = diag(Lam);
% D = diag(sum(A,2));
% L = D - A;
% [U Lam] = eig(L);
% Lam = diag(Lam);

% Crovell's filters haar wavelet
coeff1 = [1/2 -1/2];
coeff3 = [1/4 1/4 -1/4 -1/4];
coeff5 = [1/6 1/6 1/6 -1/6 -1/6 -1/6];
A0 = eye(N);

A1 = A;
D1 = diag(sum(A1,2));

Atemp = A + A^2;
A2 = double(Atemp >0);
A2 = A2-A1; % A2 is a ring adjacency matrix
D2 = diag(sum(A2,2));

Atemp = Atemp + A^3;
A3 = double(Atemp>0);
A3 = A3 - A2 - A1;
D3 = diag(sum(A3,2));

% Atemp = Atemp + A^4;
% A4 = double(Atemp>0);
% A4 = A4 - A3 - A2 - A1;
% D4 = diag(sum(A4,2));
% 
% A5 = Atemp + A^5;
% A5 = double(A5>0);
% A5 = A5 - A4 - A3 - A2 - A1;
% D5 = diag(sum(A5,2));
% 
% T1 = coeff1(1)*A0 + coeff1(2)*D1^(-1)*A1;
% T3 = coeff3(1)*A0 + coeff3(2)*D1^(-1)*A1 + coeff3(3)*D2^(-1)*A2 + coeff3(4)*D3^(-1)*A3;
% % T5 = coeff5(1)*A0 + coeff5(2)*D1^(-1)*A1 + coeff5(3)*D2^(-1)*A2 + coeff5(4)*D3^(-1)*A3 + coeff5(5)*D4^(-1)*A4 + coeff5(6)*D5^(-1)*A5;
% 
% H1 = U'*T1*U;
% H3 = U'*T3*U;
% H5 = U'*T5*U;
coeff1(1) = 10^7*quad(@max_hat,0,1/2);
coeff1(2) = 10^7*quad(@max_hat,1/2,1);

coeff3(1) = 10^3*quad(@max_hat,0,1/4);
coeff3(2) = 10^3*quad(@max_hat,0.25,0.5);
coeff3(3) = 10^3*quad(@max_hat,0.5,0.75);
coeff3(4) = 10^3*quad(@max_hat,0.75,1);

T1 = coeff1(1)*A0 + coeff1(2)*D1^(-1)*A1;
T3 = coeff3(1)*A0 + coeff3(2)*D1^(-1)*A1 + coeff3(3)*D2^(-1)*A2 + coeff3(4)*D3^(-1)*A3;
H1 = U'*T1*U;
H3 = U'*T3*U;

% crovella's filters mexican hat


% t = 0:0.01:1;
% sig = 0.25;
% phi = 2/sqrt(3*sqrt(pi)*sig)*(1-t.^2/sig^2).*exp(-t.^2/(2*sig^2));

