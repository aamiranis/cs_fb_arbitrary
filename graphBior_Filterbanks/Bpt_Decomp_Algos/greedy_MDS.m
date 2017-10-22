function [E O] = greedy_MDS(G)
% The algorithm finds an approximately minimum Dominating set in a directed
% graph using greedy heuristics
N = length(G);
% G(N,:) = 0; % do not want base staion to be selected 
% exclude the sink
min_nbr = 1; % this is the least number of even nbrs that an odd node must have
E = [];
O = [];
v_weight = ones(1,N); % assign initial weight 1 to all vertices
wnode = ones(N,1);
e_weight = (v_weight'*v_weight).*G;
outDegree = sum(e_weight,2);
inDegree = sum(e_weight);
[Min loc] = min(wnode./(outDegree+10^-10));
remain = 1:N;
% figure,

% initial check : Nodes with inDegree less than min_nbr should be assigned
% even parity
% lowDegreeNodes= find(inDegree < min_nbr);
% len = length(lowDegreeNodes);
% for i = 1:len
%     loc = lowDegreeNodes(i);
%     E = union(E,loc);
%     v_weight(loc) = 0;
%     nbr = find(e_weight(loc,:) > 0);
%     v_weight(nbr) = v_weight(nbr) - 1/min_nbr;
%     temp_index = find(v_weight(nbr) == 0);
%     O = union(O,nbr(temp_index));
%     e_weight = (v_weight'*v_weight).*G;
% end
while ~isempty(remain)
    % first check for isolated points in graph
    setA = find(v_weight == 1);
    setB = find(outDegree == 0);
    setC = find(inDegree == 0);
    isolated = intersect(intersect(setA,setB),setC);
    E = union(E, isolated);
    v_weight(isolated) = 0;
% 
%     % second check if any node partially selected has lost all its edges
%     setA = find(v_weight < 1 & v_weight > 0);
%     setB = find(outDegree == 0);
%     setC = find(inDegree == 0);
%     partial = intersect(intersect(setA,setB),setC);
%     O = union(O, partial);
%     v_weight(partial) = 0;

    % Now pick the node with min criteria and assign it to even
    E = union(E,loc);
    v_weight(loc) = 0;
    nbr = find(e_weight(loc,:) > 0);
    v_weight(nbr) = v_weight(nbr) - 1/min_nbr;
    temp_index = find(v_weight(nbr) == 0);
    O = union(O,nbr(temp_index));
    e_weight = (v_weight'*v_weight).*G;
    outDegree = sum(e_weight,2);
    inDegree = sum(e_weight);
    [Min loc] = min(wnode./(outDegree+10^-10));
    remain = setdiff(remain,[O E]);
    %         gplot(G,pos);
    %         hold on
    %         index = 1:N;
    %         x = pos(:,1);
    %         y = pos(:,2);
    %         T = text(x(index)+0.0001, y(index)+0.0001,mat2cell(index',ones(1,N),1),'color',[1 0 0]);
    %     Eloc = pos(E,:);
    %     plot(Eloc(:,1),Eloc(:,2),'s')
    %     Oloc = pos(O,:);
    %     plot(Oloc(:,1),Oloc(:,2),'o')
    %     pause
    %     hold off
end
% select the remaining ones into evens
setA = find(v_weight == 1);
setB = find(outDegree == 0);
setC = find(inDegree == 0);
isolated = intersect(intersect(setA,setB),setC);
E = union(E, isolated);
 E = setdiff(E,N);
 O = setdiff(O,N);

% % select partially selected nodes and put them in to odd set
% setA = find(v_weight < 1 & v_weight > 0);
% setB = find(outDegree == 0);
% setC = find(inDegree == 0);
% partial = intersect(intersect(setA,setB),setC);
% len = length(partial);
% O = union(O, partial);
% v_weight(partial) = 0;





