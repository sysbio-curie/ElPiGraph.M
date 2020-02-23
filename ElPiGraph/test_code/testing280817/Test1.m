%Test of work
%General constants

%nData = 1000000;
%dim = 2;
%nNodes = 50;
%X = rand(nData,dim);


X = load('test_data/tree23/tree23_inflated.data');
nData = size(X,1);
dim = size(X,2);
nNodes = 100;


% Randomly select nNodes nodes
ind = randsample(nData,nNodes);
NodePositions = X(ind,:);

tic;
XSquared = sum(X.^2,2);
toc;

tic;
[partition,dists] = PartitionData(X,NodePositions,100000,XSquared);
toc

tic;
[idx,dist] = knnsearch(NodePositions,X,'k',1);
toc;