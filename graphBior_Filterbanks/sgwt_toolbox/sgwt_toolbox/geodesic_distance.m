function D = geodesic_distance(A)
% computes the geodesic distance given adjacency matrix
n = length(A);
D = repmat(inf,n,n);

for k = 1:n
   f = k;
   s = 0;
   while ~isempty(f)
     D(k,f) = s;
     s = s+1;
     f = find(any(A(f,:),1) & D(k,:)==inf);
   end
 end