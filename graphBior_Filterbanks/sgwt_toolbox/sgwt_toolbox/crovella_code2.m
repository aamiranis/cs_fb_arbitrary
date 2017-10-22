Aam = 1000; % dimension of the field. 
N = 1000; % number of samples
Sample = floor(m*rand(N,2)) + 1;
iSample = Sample(:,1)*m+Sample(:,1)*m;
[tempSample I J] = unique(iSample);
loc = Sample(I,:);
N = length(loc);
% generate Graph
plot_option = 0;
G = graph_generator(loc,m/10);
G = double(G);
% gplot(G,loc);

D1 = geodesic_distance(G);
D2 = CT_distance(G);

for i = 1:N
    plot(sort(D2(i,:)))
    pause
end

% n = length(G);
% D1 = diag(sum(G,2));
% P = inv(D1)*G;
% [pi_ lam] = eigs(P,1);
% pi_ = pi_/sum(pi_);
% Q = ones(N,1)*pi_';
% PI_ = diag(pi_);
% P1 = P;
% 
% d1 = P1*inv(PI_);
% d2 = P1^2*inv(PI_);
% d3 = P1^3*inv(PI_);
% d4 = P1^4*inv(PI_);
% d5 = P1^5*inv(PI_);
% d6 = P1^6*inv(PI_);
% d7 = P1^7*inv(PI_);
% d8 = P1^8*inv(PI_);
% d9 = P1^9*inv(PI_);
% d10 = P1^10*inv(PI_);
% d11 = P1^11*inv(PI_);
% d12 = P1^12*inv(PI_);
% d13 = P1^13*inv(PI_);
% d14 = P1^14*inv(PI_);
% d15 = P1^15*inv(PI_);
% d16 = P1^16*inv(PI_);
% 
% for i = 1:N-1
%     for j = i+1:N
%     d = [d1(i,j) d2(i,j) d3(i,j) d4(i,j) d5(i,j) d6(i,j) d7(i,j) d8(i,j) d9(i,j)  d10(i,j)  d11(i,j)  d12(i,j)  d13(i,j)  d14(i,j)  d15(i,j)  d16(i,j) ];
%     plot(d,'*');
%     pause
%     end
% end