function I = synthetic_image()
I1 = zeros(256);
I1(1: 128,:) = 1;
I2 = triu(zeros(256)) + tril(ones(256));
% I2 = (triu(zeros(256)) + tril(ones(256)));
I2 = I2.*rot90(I2,1);
I = 0*ones(300,256);
I(1:128,:) = I1(1:128,:);
I(173:300,:) = I2(129:256,:);
I = I(100:256,85:235);
% I = ones(256);