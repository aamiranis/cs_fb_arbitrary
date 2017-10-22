max_iter = 50;
nullspace_rank = zeros(max_iter,1);
L_H = zeros(max_iter,1);
for iter = 1:max_iter
    iter
    [A xy] = random_bipartite_graph();
    F = DSATUR(A); 
    [beta bptG beta_dist Colorednodes]= harary_decomp(A,F);
    theta = size(beta,2);
    for i = 1:theta
        d1 = sum(bptG(:,:,i),2);
        d1(d1 == 0) = 1; % for isolated nodes
        d1_inv = d1.^(-0.5);
        D1_inv = diag(d1_inv);
        An = D1_inv*bptG(:,:,i)*D1_inv;
        An = 0.5*(An +An');
    end
%     P1 = nulbasis(An);
    [R, pivcol] = rref(An, sqrt(eps));
    N = length(A);
    r = length(pivcol);
    nullspace_rank(iter) = N - r;
    L_H(iter) = abs(length(Colorednodes{1}) - length(Colorednodes{2})); 
end

    
      