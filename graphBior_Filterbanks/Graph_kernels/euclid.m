function [x,y] = euclid(n,m)

% This function do the Euclid's algorithm. As a matter of fact, for two 
% given polynomials n, m (which are the polynomials of the symbolic
% variable "s") it gives two other polynomials x, y such that nx+my=1. 
% For example, use the following rules in workspace:
%
% >> syms s
% >> n = s^2;
% >> m = 6*s^2-5*s+1;
% >> [x,y]=euclid(n,m)
%
% which gives:
%
% x =
%  
% -30*s+19
%  
% y =
%  
% 1+5*s
%
% IMPORTANT NOTE: arrange n and m such that deg(n)>=deg(m).
%

syms s
r(1) = n;
r(2) = m;
q(1) = n;
x(1) = s; x(2) = 0;
x(1) = 1;
y(1) = s;
y(1) = 0; y(2) = 1;
k = 2;
while quorem(r(k),s)~=0
    k = k+1;
    [q(k-1),r(k)] = quorem(r(k-2),r(k-1));
    x(k) = x(k-2) - x(k-1)*q(k-1);
    y(k) = y(k-2) - y(k-1)*q(k-1);
end
x = x(k);
y = y(k);
% if quorem(r1,s)==0
%     x = 1/r1;
%     y = -q1/r1;
%     return
% else
%     [q2,r2] = quorem(m,r1);
%     if quorem(r2,s)==0
%         x = -q2/r2;
%         y = (1+q1*q2)/r2;
%         return
%     else
%         [q3,r3] = quorem(r1,r2);
%         x = (1+q2*q3)/r3;
%         y = (-q3-q1*(1+q2*q3))/r3;
%     end
% end
