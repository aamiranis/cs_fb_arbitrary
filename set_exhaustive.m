function [ S ] = set_exhaustive( G, H )
%SET_GREEDY Summary of this function goes here
%   Detailed explanation goes here

    n = size(H,1);
    
    S_list = nchoosek(1:n,n/2);
    min_obj = inf;
    for i = 1:size(S_list,1)
        S_curr = S_list(i,:);
        T = G(S_curr,:)'*H(S_curr,:);
        obj = norm(T - eye(size(T)),'fro');
        if (obj < min_obj)
            S = false(n,1);
            S(S_curr) = 1;
            min_obj = obj;
        end
    end
    
end
