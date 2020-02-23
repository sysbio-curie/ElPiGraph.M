display('Comparing global and local optimization on a relatively large graph and dataset');

data = load('test_data/tree23/tree23_inflated.data');

TrimmingRadius = +Inf;
NumberOfNodes = 200;

display(sprintf('Number of data points = %i',size(data,1)));

max_candidate_map = containers.Map();
max_candidate_map('bisectedge') = 10;
max_candidate_map('addnode2node') = 20;

display('Global optimization:');
tic; 
[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,NumberOfNodes,...
    'Plots',2,'verbose',1,'BranchingControls',[0.01,1],...
    'ShrinkGrammars',[],'GrowGrammars',[{'bisectedge';'addnode2node'}],...
    'TrimmingRadius',TrimmingRadius,...
    'MaxNumberOfCandidateGraphTopologiesMap',max_candidate_map); 
toc;

ln = ones(1,max(max(Edges)));
le = ones(1,size(Edges,1));
display(sprintf('Barcode for the tree = %s',getPrimitiveGraphStructureBarCode(Encode2ElasticMatrix(Edges,le,ln))));
display('Local optimization:');
%figure;
%tic; [NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,100,'Plots',2,'LocalSearch',2,'verbose',0,'BranchingControls',[0.01,1]); toc;
ln = ones(1,max(max(Edges)));
le = ones(1,size(Edges,1));
display(sprintf('Barcode for the tree = %s',getPrimitiveGraphStructureBarCode(Encode2ElasticMatrix(Edges,le,ln))));
