function y = symconv2(x, h, direction)
% symmetrically extended convolution(see section 6.5.2 in [1]):
%    x[n], E<=n<=F-1, is extended to x~[n] = x[n], 0<=n<=N-1;
%                                  x~[E-i] = x~[E+i], for all integer i
%                                x~[F-1-i] = x~[F-1+i], for all integer i
%    For odd-length h[n], to convolve x[n] and h[n], we just need extend x 
%    by (length(h)-1)/2  for both left and right edges. 
% The symmetric extension handled here is not the same as in Matlab 
%  wavelets toolbox nor in [2]. The last two use the following method:
%    x[n], E<=n<=F-1, is extended to x~[n] = x[n], 0<=n<=N-1;
%                                  x~[E-i] = x~[E+i-1], for all integer i
%                                x~[F-1-i] = x~[F+i], for all integer i 

l = length(h); s = size(x);
lext = (l-1)/2; % length of h is odd 
h = h(:)'; % make sure h is row vector 
y = x;
if strcmp(direction, 'row') % convolving along rows
    if ~isempty(x) && s(2) > 1 % unit length array skip convolution stage
        for i = 1: lext
            x = [x(:, 2*i), x, x(:, s(2)-1)]; % symmetric extension
        end
        x = conv2(x, h);
        y = x(:, l:s(2)+l-1); 
    end
elseif strcmp(direction, 'col') % convolving along columns
    if ~isempty(x) && s(1) > 1 % unit length array skip convolution stage
        for i = 1: lext 
            x = [x(2*i, :); x; x(s(1)-1, :)]; % symmetric extension
        end
        x = conv2(x, h');
        y = x(l:s(1)+l-1, :);
    end
end    
% EOF