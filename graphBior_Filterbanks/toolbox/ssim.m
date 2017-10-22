function outp = ssim(img1, img2)
f1 = uint8(img1(:));
f2 = uint8(img2(:));
f1 = double(f1);
f2 = double(f2);
mu1 = mean(f1);
mu2 = mean(f2);
var1 = var(f1);
var2 = var(f2);
var12 = cov(f1,f2);
var12 = var12(1,2);
k1 = 0.01;
k2 = 0.03;
L  = 2^8 -1;
c1 = (k1*L)^2;
c2 = (k2*L)^2;
outp = (2*mu1*mu2 +c1)*(2*var12 + c2)/((mu1^2+mu2^2 + c1)*(var1 + var2 + c2));

