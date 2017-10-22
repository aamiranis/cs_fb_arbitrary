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
%   This file implements an Meyer low-pass graph spectral kernel and
%   computes its value on the lambda values x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = meyer_kernel(x)
x = pi*x;
n = length(x);
h = zeros(n,1);
for i = 1:n
    if x(i)<=0
        h(i) = sqrt(2*theta_fun(2+1.5*x(i)/pi));
    else
        h(i) = sqrt(2*theta_fun(2-1.5*x(i)/pi));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    function theta = theta_fun(x)
        if x<=0
            theta = 0;
        elseif x>=1
            theta = 1;
        else
            theta= 3*x^2 - 2*x^3;
        end
    end
h = h';
end
