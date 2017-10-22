function y = max_hat(x,dmax)
x = x/dmax;
sigma = 0.25;
y = 2/sqrt(3*sqrt(pi)*sigma)*(1-x.^2/sigma^2).*exp(-x.^2/(2*sigma^2));
