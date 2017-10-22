function [beta bptG beta_dist Colorednodes] = twoD_decomp1(A)
% A1 = A;
theta =2;
N = length(A);
beta = zeros(N,theta);
bptG = zeros(N,N,theta);
[U Lam] = eigs(A,1,'SA');
% IDX = kmeans(U,2);
gamma = sign(U);
gamma(gamma ==0) = 1;

A_g = 1/2*(A - diag(gamma)*A*diag(gamma));
b_temp = rand(N,1);
b_temp = double(b_temp >0.5);
b_temp = 2*b_temp -1;
sats_old = 0;
B_temp = b_temp*b_temp';
c_1 = 1/2*(A - A.*B_temp)*ones(N,1);
c_2 = 1/2*(A_g - A_g.*B_temp)*ones(N,1);
sats_new_max = length(find(c_1)) + length(find(c_2));
while((sats_new_max - sats_old)>0)
    sats_old = sats_new_max;
    sats_new = zeros(N,1);
    for i = 1:N
        b_temp_temp = b_temp;
        b_temp_temp(i) = -b_temp_temp(i);
        B_temp = b_temp_temp*b_temp_temp';
        c_1 = 1/2*(A - A.*B_temp)*ones(N,1);
        c_2 = 1/2*(A_g - A_g.*B_temp)*ones(N,1);
        sats_new(i) = length(find(c_1)) + length(find(c_2));
    end
    [sats_new_max loc] = max(sats_new);
    b_temp(loc) = -b_temp(loc);
end
b_temp(loc)= -b_temp(loc);
beta(:,1) = b_temp;
beta(:,2) = gamma.*beta(:,1);
for i = 1:theta
    S1 = find(beta(:,i) == 1);
    S2 = find(beta(:,i) == -1);
    bptG(S1,S2,i) = A(S1,S2);
    bptG(S2,S1,i) = A(S2,S1);
    A(S1,S2) = 0;
    A(S2,S1) = 0;
end
Colorednodes{1} = find(beta(:,1) >0 & beta(:,2) >0);
Colorednodes{2} = find(beta(:,1) >0 & beta(:,2) <0);
Colorednodes{3} = find(beta(:,1) <0 & beta(:,2) >0);
Colorednodes{4} = find(beta(:,1) <0 & beta(:,2) <0);
beta_dist = [1 1; 1 -1; -1 1; -1 -1];

% d = sum(A1,2);
% D1 =diag(d);
% A1 = D1-A1;
% mate = card_match(A);
% S = find(mate);
% mate = [S',mate(S)'];
% 
% [U Lam] = eigs(A1,1,'SA');
% IDX1 = sign(U);
% IDX1(IDX1 ==0) = 1;
% % for i = 1:length(mate)
% %     if IDX1(mate(i,1))*IDX1(mate(i,2)) > 0
% %         [~ ,ind1] = max([U(mate(i,1)) U(mate(i,2))]);
% %         [~ ,ind2] = min([U(mate(i,1)) U(mate(i,2))]);
% %         IDX1(mate(i,ind1)) = 1;
% %         IDX1(mate(i,ind2)) = -1;
% %     end
% % end
% %         
%     
% A2 = A(Color23,Color23);
% mate = card_match(A2);
% S = find(mate);
% mate = [S',mate(S)'];
% [U Lam] = eigs(A2,1,'SA');
% IDX2 = sign(U);
% IDX2(IDX2 ==0) = 1;
% % for i = 1:length(mate)
% %     if IDX2(mate(i,1))*IDX2(mate(i,2)) > 0
% %         [~ ,ind1] = max([U(mate(i,1)) U(mate(i,2))]);
% %         [~ ,ind2] = min([U(mate(i,1)) U(mate(i,2))]);
% %         IDX2(mate(i,ind1)) = 1;
% %         IDX2(mate(i,ind2)) = -1;
% %     end
% % end
% 
% 
% % d = sum(A2,2);
% % D2 = diag(d);
% % A2 = D2 - A2;
% 
% C1 = Color14(IDX1 ==1);
% C4 = Color14(IDX1 ==-1);
% if length(C4) > length(C1)
%     C = C1;
%     C1 = C4;
%     C4 = C;
% end
% C2 = Color23(IDX2 ==1);
% C3 = Color23(IDX2 ==-1);
% % if length(C3) > length(C2)
% %     C = C2;
% %     C2 = C3;
% %     C3 = C;
% % end
% E13 = sum(sum(A(C1,C3)));
% E24 = sum(sum(A(C2,C4)));
% E12 = sum(sum(A(C1,C2)));
% E34 = sum(sum(A(C4,C3)));
% 
% if((E13 + E24) <  (E12 + E34))
%     C = C2;
%     C2 = C3;
%     C3 = C;
% end
%     
% Colorednodes{1} = C1;
% Colorednodes{2} = C2;
% Colorednodes{3} = C3;
% Colorednodes{4} = C4;
% beta_dist = [1 1; 1 -1; -1 1; -1 -1];
% beta(C1,1) = 1;
% beta(C2,1) = 1;
% beta(C3,1) = -1;
% beta(C4,1) = -1;
% beta(C1,2) = 1;
% beta(C2,2) = -1;
% beta(C3,2) = 1;
% beta(C4,2) = -1;

% 
% 
% 
% 
% % Fmin = min(F);
% % F = F - Fmin +1; % minimum color is 1
% % Fmax = max(F); % maximum value of the color
% % S = sprintf('%s%d%s','The graph has proper ', Fmax, ' coloring');
% % disp(S);
% % theta = ceil(log2(Fmax)); % this many bpts will be produced
% % disp('Computing bipartite subgraph decomposition...');
% % new_order =  dec2bin(0:(2^(theta)-1));
% % new_order = fliplr(new_order);
% % new_order = bin2dec(new_order);
% % new_order = new_order(new_order < Fmax) + 1;
% % Colorednodes = cell(Fmax,1);
% % for i = 1:Fmax
% %     Colorednodes{i} = find(F == new_order(i)); % these nodes have ith color
% % end
% % 
% % beta_dist = dec2bin(sort_colors(Fmax,theta));
% % beta_dist = double(beta_dist) - double('0');
% % beta_dist = 2*beta_dist - 1;
% % beta = zeros(N,theta);
% % bptG = zeros(N,N,theta);
% % edges = sum(sum(bptG,1),2)/2;
% % edges = edges(:);
% % cum_edges = cumsum(edges)/sum(edges);
% % theta1 = find(cum_edges > 0.97, 1,'first');
% % theta = theta1;
% % % theta =1;
% % bptG = bptG(:,:,1:theta);
% % beta = beta(:,1:theta);
% % beta_dist = beta_dist(:,1:theta);
% % beta_dist1 = 0.5*(beta_dist + 1);
% % b = 2.^((theta-1):-1:0);
% % b = b';
% % Fmax1 = 2^theta;
% % F = Fmax1 - beta_dist1*b;
% % Colorednodes1 = cell(Fmax1,1);
% % for i = 1: Fmax1
% %     Colorednodes1{i} = [];
% % end
% % for i = 1: Fmax
% %     Colorednodes1{F(i)} = union(Colorednodes1{F(i)}, Colorednodes{i});
% % end
% % Colorednodes = Colorednodes1;
% % beta_dist = dec2bin(sort_colors(Fmax1,theta));
% % beta_dist = double(beta_dist) - double('0');
% % beta_dist = 2*beta_dist - 1;
% %     
% % 
