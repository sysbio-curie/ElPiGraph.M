addpath C:\MyPrograms\_github\LizardBrain\

RandomSeed = 23;
%RandomSeed = 1;

[data,irx,v,u,s,np,ed] = lizard_brain(5,5,'make_pca_plot',0,...
    'min_branch_points',10,'epsilon',0.01,...
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
NumberOfNodes = 150;

ThresholdDistance = 0.05;

 [NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,50,...
    'Plots',2,'verbose',0,'BranchingControls',[Alpha,1],...
    'ShrinkGrammars',[],'GrowGrammars',[{'bisectedge';'addnode2node'}],...
    'TrimmingRadius',TrimmingRadius,'Lambda',Lambda,'Mu',Mu);

figure;
[dist_graph,dists,reciprocal_dists] = distance_between_nodes(np,ed,NodePositions,Edges,1,1);

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

TruePositives = length(reciprocal_dists);
FalsePositives =  NumStars_constructed-TruePositives;
FalseNegatives =  NumStars_groundtruth-TruePositives;
Precision = TruePositives/(TruePositives+FalsePositives);
Recall = TruePositives/(TruePositives+FalseNegatives);

fprintf('Identified stars percentage = %2.2f\n',100*DetectedStars/NumStars_groundtruth);
fprintf('Recall = %f\n',Recall);
fprintf('Precision = %f\n',Precision);






