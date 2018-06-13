NumberOfComputations = 10;

clear AllNodePositions;
clear AllAdjacencyMatrices;

filename = 'C:\Datas\Peshkin\demo\day22\bootstrap\bs';
NumberOfNodes = 250;
SamplingPercentage = 0.1;
TrimmingRadius = 200;

N = size(data,1);
m = size(data,2);
K = ceil(N*(1-SamplingPercentage));

AllNodePositions = zeros(NumberOfNodes,m,NumberOfComputations);
AllAdjacencyMatrices = zeros(NumberOfNodes,NumberOfNodes,NumberOfComputations);
setset

for i=1:NumberOfComputations

    sampling = randperm(N,K);
    factor = 1+(rand()-1)*0.4;
    tr = TrimmingRadius*factor;
    tic; [np,ed] = computeElasticPrincipalGraph(data(sampling,:),NumberOfNodes,'TrimmingRadius',tr,'Lambda',0.001,'Mu',0.005,'Plots',2); toc;
    drawnow;
    saveas(gcf,sprintf('%s%2i.png',filename,i));
    close all;
    
    AllNodePositions(:,:,i) = np(:,:);
    elm = MakeUniformElasticMatrix(ed,1,0);
    AllAdjacencyMatrices(:,:,i) = elm(:,:);

end

