folder = 'C:\ElPiGraph\tests\';

listing = dir(folder);

fprintf('RANDOMSEED\tTIME_FAST\tTIME_FULL\tSPEEDUP\tBAR_CODE_FAST\tBAR_CODE_FULL\tSAME_BARCODE\tMSE_FAST\tMSE_FULL\tELENERGY_FAST\tELENERGY_FULL\n');
for i=1:length(listing)
    name = listing(i).name;
    if endsWith(name,'_full.mat')
        %fprintf('%s\n',name(1:end-9));
        prefix = name(1:end-9);
        mat_full = [folder prefix '_full.mat'];
        mat_fast = [folder prefix '_fast.mat'];
        load(mat_fast);
        Elapsedtime_fast = Elapsedtime;
        MSE_fast = MSE;
        ElasticEnergy_fast = ElasticEnergy;
        NodePositions_fast = NodePositions;
        Edges_fast = Edges;
        nNodes = size(NodePositions,1);
        nEdges = size(Edges,1);
        ElasticMatrix = Encode2ElasticMatrix(Edges,ones(nEdges,1),0*ones(nNodes,1));
        BarCode_fast = getPrimitiveGraphStructureBarCode(ElasticMatrix);
        
        
        load(mat_full);
        nNodes = size(NodePositions,1);
        nEdges = size(Edges,1);
        ElasticMatrix = Encode2ElasticMatrix(Edges,ones(nEdges,1),0*ones(nNodes,1));
        BarCode = getPrimitiveGraphStructureBarCode(ElasticMatrix);

        randomseed = RandomSeed;
        fprintf('%i\t%f\t%f\t%f\t%s\t%s\t%i\t%f\t%f\t%f\t%f\n',...
            randomseed, Elapsedtime_fast,Elapsedtime,Elapsedtime/Elapsedtime_fast,...
            BarCode_fast,BarCode,strcmp(BarCode_fast,BarCode),...
            MSE_fast,MSE,ElasticEnergy_fast,ElasticEnergy);
    end
end