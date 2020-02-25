addpath C:\MyPrograms\_github\LizardBrain\

NumberOfRepeats = 100;
RandomSeeds = 1:NumberOfRepeats;

for count=1:NumberOfRepeats

fprintf('Repeat:%i\n',count);
close all;
    
RandomSeed = RandomSeeds(count);
MaxNumberOfCandidates = 20;

[data,irx,v,u,s] = lizard_brain(10,5,'make_pca_plot',0,...
    'min_branch_points',100,'epsilon',0.005,...
    'random_seed',RandomSeed);
fprintf('%i data points in %iD\n',size(data,1),size(data,2));

TrimmingRadius = +Inf;
Lambda = 0.01;
Mu = 0.1;
Alpha = 0.01;
NumberOfNodes = 100;

% =================================
% FAST MODE
% =================================

max_candidate_map = containers.Map();
max_candidate_map('bisectedge') = MaxNumberOfCandidates;
max_candidate_map('addnode2node') = MaxNumberOfCandidates;

tic; 
[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,NumberOfNodes,...
    'Plots',2,'verbose',0,'BranchingControls',[Alpha,1],...
    'ShrinkGrammars',[],'GrowGrammars',[{'bisectedge';'addnode2node'}],...
    'TrimmingRadius',TrimmingRadius,'Lambda',Lambda,'Mu',Mu,...
    'MaxNumberOfCandidateGraphTopologiesMap',max_candidate_map); 
Elapsedtime = toc;
fprintf('Time spent: %f sec\n',Elapsedtime);

nNodes = size(NodePositions,1);
nEdges = size(Edges,1);
ElasticMatrix = Encode2ElasticMatrix(Edges,ones(nEdges,1),0*ones(nNodes,1));
BarCode = getPrimitiveGraphStructureBarCode(ElasticMatrix);
fprintf('Barcode = %s\n',BarCode);

datastruct.npoints = size(data,1);
datastruct.dim = size(data,2);
datastruct.X = data;
datastruct.Weights = ones(datastruct.npoints,1);
datastruct.SquaredX = sum(data.^2,2);
datastruct.XW = datastruct.X.*datastruct.Weights;

graph.nNodes = size(NodePositions,1);
graph.NodePositions = NodePositions;
graph.Lambda = Lambda;
graph.Mu = Mu;
graph.Lambdas = Lambda*ElasticMatrix;
graph.Mus = Mu*(sum(graph.Lambdas>0)>1)';
graph.TrimmingRadius = TrimmingRadius;
graph.LocalSearch = 0;
graph.eps = 0.01;
graph.MaxMemorySize = 10000000;
graph.MaxBlockSize = graph.MaxMemorySize / graph.nNodes;

[partition, dists] = ...
    PartitionData(data, NodePositions, graph.MaxBlockSize, datastruct.SquaredX, TrimmingRadius);

partstruct.partition = partition;
partstruct.dists = dists;

[ElasticEnergy, MSE, EP, RP] = ...
    ComputePrimitiveGraphElasticEnergy(datastruct, graph, partstruct);
fprintf('ElasticEnergy = %f, MSE=%f, EP=%f, RP=%f\n',ElasticEnergy, MSE, EP, RP);

filename = sprintf('%s%i%s','C:/ElPiGraph/tests/test',count,'_fast.mat');
save(filename,'Elapsedtime','NodePositions','Edges','ReportTable','ElasticEnergy', 'MSE', 'EP', 'RP');


% =================================
% EXACT MODE
% =================================

max_candidate_map('bisectedge') = +Inf;
max_candidate_map('addnode2node') = +Inf;

tic; 
[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,NumberOfNodes,...
    'Plots',2,'verbose',0,'BranchingControls',[Alpha,1],...
    'ShrinkGrammars',[],'GrowGrammars',[{'bisectedge';'addnode2node'}],...
    'TrimmingRadius',TrimmingRadius,'Lambda',Lambda,'Mu',Mu,...
    'MaxNumberOfCandidateGraphTopologiesMap',max_candidate_map); 
Elapsedtime = toc;
fprintf('Time spent: %f sec\n',Elapsedtime);

nNodes = size(NodePositions,1);
nEdges = size(Edges,1);
ElasticMatrix = Encode2ElasticMatrix(Edges,ones(nEdges,1),0*ones(nNodes,1));
BarCode = getPrimitiveGraphStructureBarCode(ElasticMatrix);
fprintf('Barcode = %s\n',BarCode);

datastruct.npoints = size(data,1);
datastruct.dim = size(data,2);
datastruct.X = data;
datastruct.Weights = ones(datastruct.npoints,1);
datastruct.SquaredX = sum(data.^2,2);
datastruct.XW = datastruct.X.*datastruct.Weights;

graph.nNodes = size(NodePositions,1);
graph.NodePositions = NodePositions;
graph.Lambda = Lambda;
graph.Mu = Mu;
graph.Lambdas = Lambda*ElasticMatrix;
graph.Mus = Mu*(sum(graph.Lambdas>0)>1)';
graph.TrimmingRadius = TrimmingRadius;
graph.LocalSearch = 0;
graph.eps = 0.01;
graph.MaxMemorySize = 10000000;
graph.MaxBlockSize = graph.MaxMemorySize / graph.nNodes;

[partition, dists] = ...
    PartitionData(data, NodePositions, graph.MaxBlockSize, datastruct.SquaredX, TrimmingRadius);

partstruct.partition = partition;
partstruct.dists = dists;

[ElasticEnergy, MSE, EP, RP] = ...
    ComputePrimitiveGraphElasticEnergy(datastruct, graph, partstruct);
fprintf('ElasticEnergy = %f, MSE=%f, EP=%f, RP=%f\n',ElasticEnergy, MSE, EP, RP);

filename = sprintf('%s%i%s','C:/ElPiGraph/tests/test',count,'_full.mat');
save(filename,'Elapsedtime','RandomSeed','NodePositions','Edges','ReportTable','ElasticEnergy', 'MSE', 'EP', 'RP');

end