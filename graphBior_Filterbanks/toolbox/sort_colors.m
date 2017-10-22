function beta_dist = sort_colors(Fmax, theta)
if theta == 1
    beta_dist = (Fmax-1):-1:0;
    beta_dist = beta_dist(:);
else
    N1 = ceil(Fmax/2);
    N2 = floor(Fmax/2);
    beta1 = sort_colors(N1, theta -1);
    beta2 = sort_colors(N2, theta -1);
    beta_dist = [beta1(:)+2^(theta-1);beta2(:)];
end


        