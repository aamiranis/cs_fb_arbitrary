% sgwt_meshmat : Adjacency matrix for regular 2d mesh
%
% function A=meshmat_p(dim,varargin)
%
% Inputs:
% dim - size of 2d mesh
% Selectable control parameters:
% boundary - 'rectangle' or 'torus'
%
% Outputs:
% A - adjacency matrix

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond.
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

function A=sgwt_meshmat(dim,varargin)
control_params={'boundary','rectangle'};
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
if (numel(dim)==1)
    dim=[1 1]*dim;
end
% build adjacency matrix : find i,j coordinates of center points
% and right and bottom neighbors, then build connectivity matrix.
% For each valid center,neighbor pair, will add A(center,neighbor)=1
% and A(neighbor,center)=1, so A will be symmetric
N=prod(dim);
m = dim(1); n = dim(2);
[alli,allj]=find(ones(dim));
switch boundary
    case 'rectangle'
        % (ci(k),cj(k)) has neighbor (ni(k),nj(k))
        ci=[alli;alli];
        cj=[allj;allj];
        ni=[alli  ; alli+1];
        nj=[allj+1; allj];
        % prune edges at boundary
        valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2));
        ni=ni(valid);
        nj=nj(valid);
        ci=ci(valid);
        cj=cj(valid);
        cind=dim(1)*(cj-1)+ci;
        nind=dim(1)*(nj-1)+ni;
    case 'torus'
        % (ci(k),cj(k)) has neighbor (ni(k),nj(k))
        ci=[alli;alli];
        cj=[allj;allj];
        ni=[alli  ; alli+1];
        nj=[allj+1; allj];
        % wrap indices to make torus
        ni=mod(ni,dim(1))+1;
        nj=mod(nj,dim(2))+1;
        cind=dim(1)*(cj-1)+ci;
        nind=dim(1)*(nj-1)+ni;
    case 'diamond'
        % (ci(k),cj(k)) has neighbor (ni(k),nj(k))
        ci=[alli;alli];
        cj=[allj;allj];
        ni=[alli+1  ; alli+1];
        nj=[allj-1; allj+1];
        % prune edges at boundary
        valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2));
        ni=ni(valid);
        nj=nj(valid);
        ci=ci(valid);
        cj=cj(valid);
        cind=dim(1)*(cj-1)+ci;
        nind=dim(1)*(nj-1)+ni;
    case 'horizontal'
        ci=[alli];%alli];
        cj=[allj];%allj];
        ni=[alli];%  ; alli+1];
        nj=[allj+1];%; allj];
        % prune edges at boundary
        valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2));
        ni=ni(valid);
        nj=nj(valid);
        ci=ci(valid);
        cj=cj(valid);
        cind=dim(1)*(cj-1)+ci;
        nind=dim(1)*(nj-1)+ni;
    case 'diagonal'
        ci=[alli];
        cj=[allj];
        ni=[alli+1];
        nj=[allj+1];
        % prune edges at boundary
        valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2));
        ni=ni(valid);
        nj=nj(valid);
        ci=ci(valid);
        cj=cj(valid);
        cind=dim(1)*(cj-1)+ci;
        nind=dim(1)*(nj-1)+ni;
     otherwise
        error('unknown boundary option');
end
% assemble connection matrix
A=sparse([cind,nind],[nind,cind],ones(1,2*numel(ni)),N,N);

% switch boundary
%     case 'rectangle'
%         D = ones(N,1);
%         for i = 1:m
%             A(i,m+i) = 2;
%             D(i) = sqrt(2);
%             A((n-1)*m+i, (n-2)*m+i) = 2;
%             D((n-1)*m+i) = sqrt(2);
%         end
%         for j = 1:n
%             A((j-1)*m+1,(j-1)*m+2) = 2;
%             D((j-1)*m+1) = sqrt(2);
%             A(j*m,j*m-1) = 2;
%             D(j*m) = sqrt(2);
%         end
%         D(1) = 2;
%         D(m) = 2;
%         D((n-1)*m+1) = 2;
%         D(N) = 2;
%         D = spdiags(D,0,N,N);
%         A = D^(-1)*A*D;
%     case 'diamond'
%         D = ones(N,1);
%         for i = 2:m-1
%             A(i,m+i-1) = 2;
%             A(i,m+i+1) = 2;
%             D(i) = sqrt(2);
%             A((n-1)*m +i , (n-2)*m+i-1) = 2;
%             A((n-1)*m +i , (n-2)*m+i+1) = 2;
%             D((n-1)*m +i) = sqrt(2);
%         end
%         for j = 2:n-1
%             A((j-1)*m + 1, (j-2)*m +2) = 2;
%             A((j-1)*m + 1, j*m +2) = 2;
%             D((j-1)*m + 1) = sqrt(2);
%             A(j*m,(j-1)*m -1) = 2;
%             A(j*m,(j+1)*m -1) = 2;
%             D(j*m) = sqrt(2);
%         end
%         %corner points
%         A(1,m+2) = 4;
%         A(m,2*m-1) = 4;
%         A((n-1)*m+1, (n-2)*m+2) = 4;
%         A(n*m, (n-1)*m -1) = 4;
%         D(1) = 2; D(m) = 2; D((n-1)*m+1) = 2; D(n*m) = 2;
%         D = spdiags(D,0,N,N);
%         A = D^(-1)*A*D;
%     case 'horizontal'
%         D = ones(N,1);
%         for i = 1:m
%             A(i,m+i) = 2;
%             A((m-1)*m + i,(m-2)*m+i) = 2;
%             D((m-1)*m + i) = sqrt(2);
%             D(i) = sqrt(2);
%         end
%         D = spdiags(D,0,N,N);
%         A = D^(-1)*A*D;
%     otherwise
%         error('unknown boundary option');
% end





