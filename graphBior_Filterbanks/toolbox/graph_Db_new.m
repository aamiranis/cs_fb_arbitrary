function [h1 h0] = graph_Db_new(k)
% this function computes graph-Dk biorthogonal filters which have 
% k and k+1 zeors at extreme eigenvalues. 
a = zeros(1, 4*k);
for i = 1:(2*k + 2)
    a(i) = nchoosek(2*k+1,i-1);
end
A = zeros(2*k);
for i = 1:2*k
    E = 2*i;
    if (E < 2*k+1)
        S = 1;
        A(i,S:E) = fliplr(a(S:E));
    else
        S = E - 2*k;
        A(i,:) = fliplr(a(S:E));
    end
end
b0 = zeros(2*k,1);
b0(1) = 1;
b = A\b0;
b = b';
b = fliplr(b);
r = roots(b);
s =-1./r;
[dontcare index] = sort(abs(imag(s)));
s = s(index);
% # of complex roots = 2*(k-1). If (k-1) is even, then we assign the real
% root to h_1 and divide complex roots equally among h_0 and h_1. If (k-1)
% is odd, then we assign the real root to h_0 and assign the remaining root
% equally among h_1 and h_0, such that h_1 has an extra pair of conjugate
% roots than h_0. This way length of h_1 is always 2*k and length of h_0 is
% 2*k-1.
index1 = 2:4:2*k-1;
index1 = sort([index1,(index1+1)]);
index0 = setdiff(2:2*k-1,index1);
if mod(k-1,2) == 0  
    index1 = [1 index1];
else
    index0 = [1 index0];
end
s1 = s(index1);
s0 = s(index0);
h1 = zeros(1,k+1);
h1(1) = 1;
for i = 1:length(s1)
    h1 = conv(h1,[s1(i) (1-s1(i))]);
end
h1 = real(h1);
h0 = [-1 2];
for i = 1:k-1
    h0 = conv(h0,[-1 2]);
end
for i = 1:length(s0)
    h0 = conv(h0,[-s0(i) (1+s0(i))]);
end
h0 = real(h0);



% 
% 
% 
% 
% % if len is
