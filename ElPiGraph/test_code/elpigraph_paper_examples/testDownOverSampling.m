mat = load('tree_300.mat');

alpha = 0.005;
lambda = 0.03;
mu = 0.3;

n_repeats = 20;

X = mat.X; 

Npoints = size(X,1)
dim = size(X,2)

maxNumNodes = 50;

figure;
[np,ed,ReportTable] = computeElasticPrincipalGraph(X,maxNumNodes,'BranchingControls',[alpha 1],'Lambda',lambda,'Mu',mu,'Plots',0); 
plot(X(:,1),X(:,2),'kx','MarkerSize',5); 
drawGraph2D(np,ed,'LineColor','r','LineWidth',3,'ShowClusterNumbers',0,'NOdeSizes',1*ones(1,size(np,1))); 
drawnow;

%figure;
clear dist_graph_down;
for j=1:5+1
    figure;
for i=1:n_repeats
subsampleSize = j*50;
if subsampleSize==300
    subsampleSize = 290;
end
dX = X(randsample(Npoints,subsampleSize),:);
[np_down,ed_down] = computeElasticPrincipalGraph(dX,maxNumNodes,'BranchingControls',[alpha 1],'Lambda',lambda,'Mu',mu,'Plots',0,'verbose',0); 
%plot(dX(:,1),dX(:,2),'k.','MarkerSize',5); 
%drawGraph2D(np_down,ed_down,'LineColor','g'); 
%drawnow;
dist_graph_down(j,i) = distance_between_nodes(np,ed,np_down,ed_down,1);
title(sprintf('Subsample size = %i',subsampleSize),'FontSize',14);
drawnow;
disp(sprintf('%i,%i: %f',size(dX,1),i,dist_graph_down(j,i)));
end
end

clear dist_graph_over;
for jj=1:10
    figure;
for ii=1:n_repeats/2

    inflationFactor = jj*2;
    X1 = zeros(size(X,1)*inflationFactor,size(X,2));
    %disp(size(X1,1));
    k=1;
    STDV = std(X);
    if 1
        for i=1:size(X,1)
            for j=1:inflationFactor
                r = (0.5-rand(1,size(X,2)))*2;
                p = 0.4*STDV.*r;
                X1(k,:) = X(i,:)+p;
                k=k+1;
            end
        end
        %X = X1;
    end
    
    [np_over,ed_over,ReportTable] = computeElasticPrincipalGraph(X1,maxNumNodes,'BranchingControls',[alpha 1],'Lambda',lambda,'Plots',0,'verbose',0); 
    dist_graph_over(jj,ii) = distance_between_nodes(np,ed,np_over,ed_over,1);
    title(sprintf('Subsample size = %i',size(X1,1)),'FontSize',14);
    drawnow;
    disp(sprintf('%i,%i: %f',size(X1,1),ii,dist_graph_over(jj,ii)));

end
end
%figure;
%[np_over,ed_over,ReportTable] = computeElasticPrincipalGraph(X,maxNumNodes,'BranchingControls',[alpha 1],'Lambda',lambda,'Plots',0); 
%plot(X(:,1),X(:,2),'k.','MarkerSize',5); 
%drawGraph2D(np_over,ed_over,'Color','r'); 
%drawnow;


figure;
n = size(dist_graph_down,1); semilogx([1:n]*50/300,mean(dist_graph_down'),'ro-','LineWidth',5); hold on; 
semilogx([1:n]*50/300,dist_graph_down','kx','LineWidth',2); 
set(gca,'FontSize',14); xlabel('Down/Over-sampling','FontSize',14); ylabel('Distance to the reference','FontSize',14); hold on;

n = size(dist_graph_over,1); semilogx([1:n]*4,mean(dist_graph_over'),'bo-','LineWidth',5); hold on; 
semilogx([1:n]*4,dist_graph_over','kx','LineWidth',2); 
set(gca,'FontSize',14); xlabel('Down/Over-sampling','FontSize',14); ylabel('Distance to the reference','FontSize',14);
