function [Ta Td] = Random_tx(A,level)
% Ramchandran transform
a = 1/4;
b =1/4;
N = length(A);
Al = A^level;
Al = Al - diag(diag(Al));
Al = double(Al>0);
Dl = diag(sum(Al,2));
Ll = Dl - Al;
Ta = eye(N) - a*(eye(N) + Dl)^(-1)*Ll;
Td = eye(N) + b*(eye(N) + Dl)^(-1)*Ll;
