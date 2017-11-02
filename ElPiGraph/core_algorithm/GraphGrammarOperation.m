function [NodePositionArray, ElasticMatrices, NodeIndicesArray] =...
    GraphGrammarOperation(NodePositions, ElasticMatrix, X, type, partition)
%
% This is the core function for application of graph grammar approach for
% constructing primitive elastic principal graphs
% 
% The function takes a definition of the principal graph embeddment and
% applies a graph grammar operation of type type
% 
% Input arguments:
%   NodePositions - position of nodes of the primitive elastic graph
%   ElasticMatrix - matrix with elasticity coefficients (lambdas for edges,
%       and star elasticities along the diagonal
%   X - is dataset which is to be approximated by the graph
%   type - one of the operation types:
%       'addnode2node'      adds a node to each graph node
%       'removenode'        removes terminal node
%       'bisectedge'        adds nodt to the middle of each edge
%       'shrinkedge'        removes edge and glues two ends of this edge
%   partition is n-by-1 vector. partition(i) is number of node which is
%       associated with data point X(i,:).
%
% Outputs:
%   NodePositionArray - 3d array with dimensions
%       [node_number, NodePosition, graph_number], represents all generated
%       node configurations.
%   ElasticMatrices - 3d array with dimensions 
%       [node_number, node_number, graph_number], represents all generated
%       elasticity matrices  
%   NodeIndicesArray in the version 1.1 with a possibility of local search,
%       each operation reports NodeIndices which specifies how the nodes in
%       the newly generated graph are related to the nodes in the initial
%       graph, and contains zeros for new nodes. 
%
    switch type
        case 'addnode2node'
            [NodePositionArray, ElasticMatrices, NodeIndicesArray] =...
                AddNode2Node(NodePositions, ElasticMatrix, X, partition);
        case 'removenode'
            [NodePositionArray, ElasticMatrices, NodeIndicesArray] =...
                RemoveNode(NodePositions, ElasticMatrix);
        case 'bisectedge'
            [NodePositionArray, ElasticMatrices, NodeIndicesArray] =...
                BisectEdge(NodePositions, ElasticMatrix);
        case 'shrinkedge'
            [NodePositionArray, ElasticMatrices, NodeIndicesArray] =...
                ShrinkEdge(NodePositions, ElasticMatrix);
        otherwise
            error('ERROR: operation %s is not defined',type);
    end
end


function [NodePositionArray, ElasticMatrices, NodeIndicesArray]...
    = AddNode2Node(NodePositions, ElasticMatrix, X, partition)
%
% This grammar operation adds a node to each graph node
% The positions of the node is chosen as a linear extrapolation for a leaf
% node (in this case the elasticity of a newborn star is chosed as in
% BisectEdge operation),
% or
% as the data point giving the minimum local MSE for a star (without any optimization).

    NNodes = size(NodePositions,1);
    NNp1 = NNodes + 1;
    NumberOfGraphs = NNodes;
    NodePositionArray = zeros(NNp1, size(NodePositions, 2), NumberOfGraphs);
    ElasticMatrices = zeros(NNp1, NNp1, NumberOfGraphs);
    NodeIndicesArray = zeros(NNp1, NumberOfGraphs);
    Mus = diag(ElasticMatrix);
    L = ElasticMatrix - diag(Mus);
    indL = L > 0;
    Connectivities = sum(indL);
    assoc = accumarray(partition + 1, 1, [NNp1, 1]);
    assoc = assoc(2:end);

    % Create prototypes for new NodePositions, ElasticMatrix and inds
    NPProt = zeros(NNp1, size(NodePositions, 2));
    NPProt(1:NNodes,:) = NodePositions;
    EMProt = zeros(NNp1, NNp1);
    EMProt(1:NNodes, 1:NNodes) = ElasticMatrix;
    NIProt = [1:NNodes,0];
    
    % Main loop
    for i=1:NNodes
        % Compute mean edge elastisity for edges with node i 
        meanLambda = mean(L(i,indL(i,:)));
        % Put prototypes to corresponding places
        NodePositionArray(:, :, i) = NPProt;
        ElasticMatrices(:, :, i) = EMProt;
        NodeIndicesArray(:, i) = NIProt;
        % Add edge to elasticity matrix
        ElasticMatrices(NNp1, i, i) = meanLambda;
        ElasticMatrices(i, NNp1, i) = meanLambda;

        if(Connectivities(i)==1)
            % Add node to terminal node
            ineighbour = find(indL(i,:));
            % Calculate new node position
            NewNodePosition = 2 * NodePositions(i, :)...
                - NodePositions(ineighbour, :);
            % Complete elasticity matrix
            ElasticMatrices(i, i, i) = Mus(ineighbour);
        else
            % Add node to star
            % If number of data points associated with star centre is zero
            if assoc(i) == 0
                % then select the mean of all leaves as new position
                NewNodePosition = mean(NodePositions(indL(:, i), :));
            else
                % otherwise take the mean of points associated with central
                % node.
                NewNodePosition = mean(X(partition == i, :));
            end 
        end
        % Fill NodePosition
        NodePositionArray(NNp1,:, i) = NewNodePosition;
    end
end


function [NodePositionArray, ElasticMatrices, NodeIndicesArray] =...
    BisectEdge(NodePositions, ElasticMatrix)
%
% This grammar operation inserts a node inside the middle of each edge
% The elasticity of the edges do not change
% The elasticity of the newborn star is chosen as 
% mean over the neighbour stars if the edge connects two star centers
% or 
% the one of the single neigbour star if this is a dangling edge
% or 
% if one starts from a single edge, the star elasticities should be on
% one of two elements in the diagoal of the ElasticMatrix 

    % Decompose Elastic Matrix: Mus
    Mus = diag(ElasticMatrix);
    % Get list of edges
    [start, stop] = find(triu(ElasticMatrix,1));
    % Define some constants
    NumberOfGraphs = length(start);
    NNodes = size(NodePositions,1);
    NNp1 = NNodes + 1;
    % Preallocate arrays
    NodePositionArray = zeros(NNp1, size(NodePositions, 2), NumberOfGraphs);
    ElasticMatrices = zeros(NNp1, NNp1, NumberOfGraphs);
    NodeIndicesArray = zeros(NNp1, NumberOfGraphs);
    % Create prototypes for new NodePositions, ElasticMatrix and inds
    NPProt = zeros(NNp1, size(NodePositions, 2));
    NPProt(1:NNodes, :) = NodePositions;
    EMProt = zeros(NNp1, NNp1);
    EMProt(1:NNodes, 1:NNodes) = ElasticMatrix;
    NIProt = [1:NNodes,0];
    % Main loop
    for i=1:NumberOfGraphs
        NewNodePosition = (NodePositions(start(i), :)...
            + NodePositions(stop(i), :)) / 2;
        % Put prototypes to corresponding places
        NodePositionArray(:, :, i) = NPProt;
        ElasticMatrices(:, :, i) = EMProt;
        NodeIndicesArray(:, i) = NIProt;
        % Fill NodePosition
        NodePositionArray(NNp1,:, i) = NewNodePosition;
        %Correct ElasticMatrix
        lambda = ElasticMatrix(start(i), stop(i));
        % Remove edge
        ElasticMatrices(start(i), stop(i), i) = 0;
        ElasticMatrices(stop(i), start(i), i) = 0;
        % Add two edges
        ElasticMatrices(start(i), NNp1, i) = lambda;
        ElasticMatrices(NNp1, start(i), i) = lambda;
        ElasticMatrices(NNp1, stop(i), i) = lambda;
        ElasticMatrices(stop(i), NNp1, i) = lambda;
        % Define Mus of edge nodes:
        mu1 = Mus(start(i));
        mu2 = Mus(stop(i));
       if mu1 > 0 && mu2 > 0
           ElasticMatrices(NNp1, NNp1, i) = (mu1 + mu2) / 2;
       else
           ElasticMatrices(NNp1, NNp1, i) = max(mu1, mu2);
       end
    end
end


function [NodePositionArray, ElasticMatrices, NodeIndicesArray] =...
    RemoveNode(NodePositions, ElasticMatrix)
%
% This grammar operation removes a leaf node (connectivity==1)
%
    % Decompose ElasticMatrix
    Mus = diag(ElasticMatrix);
    L = ElasticMatrix - diag(Mus);
    Connectivities = sum(L > 0);
    % Define sizes
    NNodes = size(NodePositions,1);
    NumberOfGraphs = sum(Connectivities==1);
    % Preallocate arrays
    NodePositionArray = zeros(NNodes-1, size(NodePositions,2), NumberOfGraphs);
    ElasticMatrices = zeros(NNodes-1, NNodes-1, NumberOfGraphs);
    NodeIndicesArray = zeros(NNodes-1, NumberOfGraphs);
    k=1;
    for i=1:length(Connectivities)
        if Connectivities(i) == 1
            % It is terminal node. Remove it
            newinds = [1:(i - 1),(i + 1):NNodes];
            NodePositionArray(:, :, k) = NodePositions(newinds, :);
            ElasticMatrices(:, :, k) = ElasticMatrix(newinds, newinds);
            NodeIndicesArray(:, k) = newinds;
            k=k+1;
        end
    end
end

function [NodePositionArray, ElasticMatrices, NodeIndicesArray] = ShrinkEdge(NodePositions,ElasticMatrix)
%
% This grammar operation removes an edge from the graph
% If this is an edge connecting a leaf node then it is equivalent to
% RemoveNode. So we remove only internal edges.
% If this is an edge connecting two stars then their leaves are merged,
% and the star is placed in the middle of the shrinked edge.
% The elasticity of the new formed star is the average of two star
% elasticities.
%
    % Decompose ElasticMatrix
    Mus = diag(ElasticMatrix);
    L = ElasticMatrix - diag(Mus);
    Connectivities = (sum(L > 0))';
    % Get list of edges
    [start, stop] = find(triu(ElasticMatrix,1));
    % Define sizes
    NNodes = size(NodePositions,1);
    %Identify edges with minimal connectivity greater than 1
    ind = min([Connectivities(start'),Connectivities(stop')],[],2);
    ind = ind > 1;
    start = start(ind);
    stop = stop(ind);
    %Calcualte number of graphs
    NumberOfGraphs = length(start);
    % Preallocate arrays
    NodePositionArray = zeros(NNodes-1, size(NodePositions,2), NumberOfGraphs);
    ElasticMatrices = zeros(NNodes-1, NNodes-1, NumberOfGraphs);
    NodeIndicesArray = zeros(NNodes-1, NumberOfGraphs);

    for i=1:length(stop)
        %Create copy of elastic matrix
        em = ElasticMatrix;
        % Reattaches all edges connected with stop(i) to start(i)
        % and make a new star with an elasticity average of two merged stars
        em(start(i),:) = max(L(start(i), :), L(stop(i), :));
        em(:,start(i)) = max(L(:, start(i)), L(:, stop(i)));
        em(start(i),start(i)) = (Mus(start(i)) + Mus(stop(i))) / 2;
        % Create copy of node positions
        np = NodePositions;
        % Modify node start(i)
        np(start(i), :) = (np(start(i), :) + np(stop(i), :)) / 2;
        % Form index for retained nodes and extract corresponding part of 
        % node positions and elastic matrix
        newinds = [1:(stop(i)-1),(stop(i) + 1):NNodes];
        NodePositionArray(:, :, i) = np(newinds, :);
        ElasticMatrices(:, :, i) = em(newinds, newinds);
        NodeIndicesArray(:, i) = newinds;
    end
end
