function qmf_filt = johnston_kernel1()
close all
% computes johnston kernel of length filterlen
% g = @(x)(meyer_kernel(x));  % Meyer Wavelet Kernel
h = cell(10,1);
for taps = 4:4:20
    m = taps/2;
    ind = taps/4;
    x0 = ones(m,1);
    if m >1
        x0 = 0*x0;
        x0(3:end) =h{ind-1}(1:2:end);
    end
    [xopt,Opt,Nav]=hooke(@obj_f,x0);
    xopt = xopt(:)';
    h{ind} = zeros(1,taps);
    h{ind}(1:m) = xopt;
    h{ind}((m+1):2*m) = fliplr(h{m}(1:m));
%     h{ind} = flipud(h{ind});
end
lam = 0:0.01:2;
for i = 1:ind
    figure,plot(lam,abs(polyval(h{i},lam)))
end
qmf_filt = h; 
end
function c = obj_f(x)
m = length(x);
lam = 1:0.01:2;
h = zeros(1,2*m);
h(1:m) = x(:)';
h((m+1):2*m) = fliplr(h(1:m));
S = sum(polyval(conv(h,h),lam));
lam1 = 0:0.01:2;
lam2 = fliplr(lam1);
E = 2 - polyval(conv(h,h),lam1) - polyval(conv(h,h),lam2);
E = norm(E)^2;
c = 0.5*S + 0.5*E;
end

