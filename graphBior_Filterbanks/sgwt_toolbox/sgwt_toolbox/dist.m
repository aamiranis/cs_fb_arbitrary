function D = dist(X)
N = length(X);
S = X*X';
D1 = diag(S);
D1 = diag(D1);
J = ones(N,N);

D = D1*J + J*D1 - 2*S;
D = D.^(0.5);

