function [ S ] = set_greedy_deterministic_disjoint( G, H )
%SET_GREEDY Summary of this function goes here
%   Detailed explanation goes here

    N = size(H,2);
    
    Sc = true(2*N,1);
    
    % This is the true S, updated separately from Sc
    S = false(2*N,1);
    
    q = zeros(2*N,1);
    for i = 1:2*N
        q(i) = (G(i,:)*G(i,:)') * (H(i,:)*H(i,:)'); 
    end
    
    p = zeros(2*N,1);
    for i = 1:2*N
        p(i) = - 2 * G(i,:)*H(i,:)';
    end
    
    for iter = 1:N
        indices = find(Sc);
        
        delta = (p(Sc) + q(Sc));
        
        [~,min_ind] = min(delta);
        u = indices(min_ind);
        
        % Add u to sampling set
        Sc(u) = 0;
        S(u) = 1;
        
        % Remove "images": make corresponding node in the other set unavailable
        if (u > N)
            Sc(u - N) = 0;
        else
            Sc(u + N) = 0;
        end
        
        p(Sc) = p(Sc) + 2 * (G(Sc,:)*G(u,:)').*(H(Sc,:)*H(u,:)');
    end
    
end
