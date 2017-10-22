function [ S ] = set_greedy_deterministic( G, H )
%SET_GREEDY Summary of this function goes here
%   Detailed explanation goes here

    N = size(H,2);
    Sc = true(2*N,1);
    
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
        
        Sc(u) = 0;
        
        p(Sc) = p(Sc) + 2 * (G(Sc,:)*G(u,:)').*(H(Sc,:)*H(u,:)');
    end
    
    S = ~Sc;
end
