function A = kregularize(B)
A = 0*B;
%corner points extension
s = length(B);
s = sqrt(s);
% A(1,2) = 1; A(2,1) = 1;
% A(1,s) = 1; A(s,1) = 1;
% 
% A(s, s-1) = 1; A(s-1, s) = 1; 
% A(s, 2*s) = 1; A(2*s, s) = 1;
% 
% A((s-1)*s+1,(s-2)*s+1) = 1; A((s-2)*s+1,(s-1)*s+1) = 1;
% A((s-1)*s+1,(s-1)*s+2) = 1; A((s-1)*s+2,(s-1)*s+1) = 1;
% 
% A(s*s,s*s-1) = 1; A(s*s -1,s*s) = 1;
% A(s*s,s*(s-1)) = 1; A(s*(s-1),s*s) = 1;


% other boundry points 

% left
for i = 1:2:s-1
    A(i,i+1) = 1; A(i+1,i) = 1;
end

% right
j = s*(s-1);
for i = 1:2:s-1
    A(j+i,j+i+1) = 1; A(j+i+1,j+i) = 1;
end
    

%top
j = 1; 
for i = 1:2:s-1
    A((i-1)*s+j,i*s+j) = 1; A(i*s+j,(i-1)*s+j) = 1;
end

%bottom
j = s; 
for i = 1:2:s-1
    A((i-1)*s+j,i*s+j) = 1; A(i*s+j,(i-1)*s+j) = 1;
end
A = A+B;