n = 30;
N = n^2;
[A D] = sgwt_meshmat(n);
A = full(A); D= full(D);
A = 0.5*(A+A');
L = 4*eye(N) - A;
[U Lam] = eig(L);
Lam = diag(Lam);
[Lam I ] = sort(Lam);
U = U(:,I);
g_low = [ones(1,N/2) zeros(1,N/2)];
g_high = fliplr(g_low);
T_low = U*diag(g_low)*U';
T_high = U*diag(g_high)*U';
alpha = [ones(1,N/2-5) zeros(1,N/2+5)]';
% f = U*alpha;
f = ones(N,1);
% f = f';
% f = D*f;
f_high = T_high*f;
f_low = T_low*f;
E = [];
O= [];
s_im = n;
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

f_high(O) = 0;
f_low(O) = 0;
f_high_r = 2*T_high*f_high;
% f_high_r = D^(-1)*f_high_r;
f_low_r = 2*T_low*f_low;
% f_low_r = D^(-1)*f_low_r;

figure,
plot(f);
hold on
plot(f_low_r,'r');
plot(f_high_r,'m');
