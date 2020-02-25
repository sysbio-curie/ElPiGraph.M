addpath C:\MyPrograms\_github\LizardBrain\

RandomSeed = 23;
%RandomSeed = 1;

[data,irx,v,u,s,np,ed] = lizard_brain(10,5,'make_pca_plot',0,...
    'min_branch_points',100,'epsilon',0.001,...
    'random_seed',RandomSeed,'make_pca_plot',1,'add_noise',0.03);
fprintf('%i data points in %iD\n',size(data,1),size(data,2));

nNodes = size(np,1);
nEdges = size(ed,1);
ElasticMatrix = Encode2ElasticMatrix(ed,ones(nEdges,1),0*ones(nNodes,1));
BarCode = getPrimitiveGraphStructureBarCode(ElasticMatrix);
fprintf('Ground truth barcode:%s\n',BarCode);


TrimmingRadius = 0.2;
Lambda = 0.001;
Mu = 0.05;
Alpha = 0.01;
NumberOfNodes = 50;
ThresholdDistance = 0.05;

 [NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,NumberOfNodes,...
    'Plots',2,'verbose',1,'BranchingControls',[Alpha,1],...
    'ShrinkGrammars',[],'GrowGrammars',[{'bisectedge';'addnode2node'}],...
    'TrimmingRadius',TrimmingRadius,'Lambda',Lambda,'Mu',Mu);

figure;
[dist_graph,dists,reciprocal_dists,TP,FP,FN] = distance_between_nodes(np,ed,NodePositions,Edges,1,1);

fprintf('Distance = %f\n',dist_graph);
dists
reciprocal_dists

        nNodes = size(NodePositions,1);
        nEdges = size(Edges,1);
        ElasticMatrix_constructed = Encode2ElasticMatrix(Edges,ones(nEdges,1),0*ones(nNodes,1));
        BarCode_constructed = getPrimitiveGraphStructureBarCode(ElasticMatrix_constructed);

fprintf('Constructed barcode:%s\n',BarCode_constructed);

NumStars_groundtruth = sum(sum(ElasticMatrix>0)>2);
NumStars_constructed = sum(sum(ElasticMatrix_constructed>0)>2);

Precision = TP/(TP+FP);
Recall = TP/(TP+FN);

fprintf('Recall = %f\n',Recall);
fprintf('Precision = %f\n',Precision);






