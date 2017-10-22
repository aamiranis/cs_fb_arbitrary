%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   "Copyright (c) 2011 The University of Southern California"
%   All rights reserved.
%
%   Permission to use, copy, modify, and distribute this software and its
%   documentation for any purpose, without fee, and without written
%   agreement is hereby granted, provided that the above copyright notice,
%   the following two paragraphs and the author appear in all copies of
%   this software.
%
%   NO REPRESENTATIONS ARE MADE ABOUT THE SUITABILITY OF THE SOFTWARE
%   FOR ANY	PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED
%   WARRANTY.
%
%   Neither the software developers, the Compression Research Group,
%   or USC, shall be liable for any damages suffered from using this
%   software.
%
%   Author: Sunil K Narang
%   Director: Prof. Antonio Ortega
%   Compression Research Group, University of Southern California
%   http://biron.usc.edu/wiki/index.php?title=CompressionGroup
%   Contact: kumarsun@usc.edu
%
%   Date last modified:	07/05/2011 kumarsun
%
%   Description:
%   This file implements an ideal low-pass graph spectral kernel and
%   computes its value on the lambda values x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = ideal_kernel(x)
n = length(x);
h = ones(n,1);
for i = 1:n;
    if x(i)<(1-10^(-12))
        h(i) = sqrt(2);
    elseif x(i)>(1 + 10^12)
        h(i) = 0;
    else
        h(i) = 1;
    end
end
h = h';
% delta = 0.05;
% gamma = 4*delta;
% theta = atan(1/delta);
% x0 = 1-gamma;
% y0 = cot(theta/2)*x0 + (2-cot(theta/2)*(1-delta));
% rad = 2-y0;
% x1= x0 + rad*sin(theta);
% 
% z0 = 1+gamma;
% w0 = cot(theta/2)*z0 -cot(theta/2)*(1+delta);
% rad1 = w0;
% z1 = z0 -rad1*sin(theta);
% 
% % x1 = 0:0.01:0.95;
% % y1 = sqrt(2)*ones(length(x1),1);
% % x2 = 1.05:0.01:2;
% % y2 = zeros(length(x2),1);
% % x_t = [x1 1 x2];
% % y_t = [y1 ; 1; y2]';
% for i = 1:n;
%     if x(i)<1-gamma
%         h(i) = 2;
%     elseif x(i)>1+gamma
%         h(i) = 0;
%     elseif x(i) < x1
%         h(i) = y0+rad*sin(acos((x(i)-x0)/rad));
%     elseif x(i) > z1
%         h(i) = w0 - rad1*sin(acos((x(i)-z0)/rad1));
%     else
%         h(i) = -(1/delta)*x(i)+(1+1/delta);
%     end
% end
% h2 = h;
% h = h2.^(0.5)';
