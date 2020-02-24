function [dist_graph,dists,reciprocal_dists] = distance_between_nodes(np1,ed1,np2,ed2,makeplots,only_between_stars)
% Computes distance between two graphs

if nargin<5
    makeplots = 0;
end
if nargin<6
    only_between_stars = 0;
end


dist_graph = 0;

adj_mat1 = zeros(size(np1,1),size(np1,1));
for i=1:size(ed1,1)
    adj_mat1(ed1(i,1),ed1(i,2))=1;
    adj_mat1(ed1(i,2),ed1(i,1))=1;
end
adj_mat2 = zeros(size(np2,2),size(np2,2));
for i=1:size(ed2,1)
    adj_mat2(ed2(i,1),ed2(i,2))=1;
    adj_mat2(ed2(i,2),ed2(i,1))=1;
end
conn1 = sum(adj_mat1);
conn2 = sum(adj_mat2);

if only_between_stars
    np1s = np1(conn1>2,:);
    np2s = np2(conn2>2,:);
    [partition1,dists1] = knnsearch(np1s,np2s,'k',1);
    [partition2,dists2] = knnsearch(np2s,np1s,'k',1);
else
    [partition1,dists1] = knnsearch(np1,np2,'k',1);
    [partition2,dists2] = knnsearch(np2,np1,'k',1);
end

reciprocal_pairs = [];
k =1;
for i=1:length(partition1)
    j = partition1(i);
    if i==partition2(j)
        reciprocal_pairs(k,1) = i;
        reciprocal_pairs(k,2) = j;
        k = k+1;
    end
end

% compute total length of the reference (first) tree

tree_length = 0;
for i=1:size(ed1,1)
    tree_length = tree_length + sqrt(sum((np1(ed1(i,1),:)-np1(ed1(i,2),:)).^2));
end
%disp(tree_length);

labels = {};
for i=1:10000 labels{i} = int2str(i); end

if makeplots

if size(np1,2)>2 
npconc = [np1;np2];
[v,u,s] = pca(npconc);
npu1 = u(1:size(np1,1),1:2);
npu2 = u(size(np1,1)+1:size(np2,1)+size(np1,1),1:2);
drawGraph2D(npu1,ed1,'LineColor','r','LineWidth',5,'ShowClusterNumbers',0); hold on;
drawGraph2D(npu2,ed2,'LineColor','b','LineWidth',1,'ShowClusterNumbers',0); hold on;
plot(npu1(conn1>2,1),npu1(conn1>2,2),'ro','MarkerSize',8,'LineWidth',5);
text(npu1(conn1>2,1)+0.05,npu1(conn1>2,2),labels(conn1>2),'Color','k');
plot(npu2(conn2>2,1),npu2(conn2>2,2),'bo','MarkerSize',5);
text(npu2(conn2>2,1),npu2(conn2>2,2),labels(conn2>2),'Color','b');
else
drawGraph2D(np1,ed1,'LineColor','r','LineWidth',5,'ShowClusterNumbers',0); hold on;
drawGraph2D(np2,ed2,'LineColor','b','LineWidth',1,'ShowClusterNumbers',0); hold on;
plot(np1(conn1>2,1),np1(conn1>2,1),'ro','MarkerSize',5);
plot(np2(conn2>2,1),np2(conn2>2,1),'bo','MarkerSize',5);
end
end

if only_between_stars

starind1 = find(conn1>2);
starind2 = find(conn2>2);
    
for j=1:length(partition1)
    j1 = j;
    i = partition1(j1);
    j1 = starind2(j1);
    i1 = starind1(i);
    if makeplots
        if size(np1,2)==2 
            plot([np1(i1,1) np2(j1,1)],[np1(i1,2) np2(j1,2)],'k--','LineWidth',0.5,'Color',[0.4 0.4 0.4]);
        else
            plot([npu1(i1,1) npu2(j1,1)],[npu1(i1,2) npu2(j1,2)],'k--','LineWidth',0.5,'Color',[0.4 0.4 0.4]);
        end
    end
    dists(j) = sqrt(sum((np1(i1,:)-np2(j1,:)).^2))/tree_length;
    dist_graph = dist_graph+sum((np1(i1,:)-np2(j1,:)).^2);
    
    %if partition2(partition1(i))==i
        %plot([np1(i,1) np2(partition1(i),1)],[np1(i,2) np2(partition1(i),2)],'k-');
    %end
end    

for k=1:size(reciprocal_pairs,1)
    i = reciprocal_pairs(k,2);
    j = reciprocal_pairs(k,1);
    i = starind1(i);
    j = starind2(j);
    if makeplots
        if size(np1,2)==2 
            plot([np1(i,1) np2(j,1)],[np1(i,2) np2(j,2)],'k-','LineWidth',2);
        else
            plot([npu1(i,1) npu2(j,1)],[npu1(i,2) npu2(j,2)],'k-','LineWidth',2);
        end
    end
    reciprocal_dists(k) = sqrt(sum((np1(i,:)-np2(j,:)).^2))/tree_length;
    
end
    
else
    
for j=1:length(partition1)
    j1 = j;
    i = partition1(j1);
    if makeplots
        if size(np1,2)==2 
            plot([np1(i,1) np2(j,1)],[np1(i,2) np2(j,2)],'k--','LineWidth',0.5,'Color',[0.4 0.4 0.4]);
        else
            plot([npu1(i,1) npu2(j,1)],[npu1(i,2) npu2(j,2)],'k--','LineWidth',0.5,'Color',[0.4 0.4 0.4]);
        end
    end
    dists(j) = sqrt(sum((np1(i,:)-np2(j,:)).^2))/tree_length;
    dist_graph = dist_graph+sum((np1(i,:)-np2(j,:)).^2);
    %if partition2(partition1(i))==i
        %plot([np1(i,1) np2(partition1(i),1)],[np1(i,2) np2(partition1(i),2)],'k-');
    %end
end
end

dist_graph = sqrt(dist_graph)/tree_length;

end

