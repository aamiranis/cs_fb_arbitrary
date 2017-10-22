function y = haar_coeffs(x,scale)
pos1 = find(x*2^scale <= 0.5);
pos2 = find(x*2^scale > 0.5);
y = x;
y(pos1) = 0.5;
y(pos2) = -0.5;

