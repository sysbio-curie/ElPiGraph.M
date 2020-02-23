function [dist_graph] = distance_between_nodes(np1,ed1,np2,ed2,makeplots)
% Computes distance between two graphs

if nargin<5
    makeplots = 0;
end

dist_graph = 0;

[partition1,dists1] = knnsearch(np1,np2,'k',1);
[partition2,dists2] = knnsearch(np2,np1,'k',1);

% compute total length of the reference (first) tree

tree_length = 0;
for i=1:size(ed1,1)
    tree_length = tree_length + sqrt((np1(ed1(i,1),1)-np1(ed1(i,2),1))^2+(np1(ed1(i,1),2)-np1(ed1(i,2),2))^2);
end
%disp(tree_length);

if makeplots
drawGraph2D(np1,ed1,'LineColor','r','LineWidth',5,'ShowClusterNumbers',0); hold on;
drawGraph2D(np2,ed2,'LineColor','b','LineWidth',1,'ShowClusterNumbers',0); hold on;
end

for j=1:length(partition1)
    i = partition1(j);
    if makeplots
    plot([np1(i,1) np2(j,1)],[np1(i,2) np2(j,2)],'k--','LineWidth',0.5,'Color',[0.4 0.4 0.4]);
    end
    dist_graph = dist_graph+((np1(i,1)-np2(j,1))^2+(np1(i,2)-np2(j,2))^2);
    %if partition2(partition1(i))==i
        %plot([np1(i,1) np2(partition1(i),1)],[np1(i,2) np2(partition1(i),2)],'k-');
    %end
end

dist_graph = sqrt(dist_graph)/tree_length;

end

