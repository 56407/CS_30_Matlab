function nodeIds = extractNodeIdsOnEdges(e, edgeIds)

% extractNodeIdsOnEdges - Extract ids of nodes lying on given edges.
%
% This QuickerSim CFD Toolbox function is intended for convenient
% preparation of plots of solution along boundaries of the domain. It
% allows to quickly find ids of the nodes which are lying on the given
% edges. Later this ids may be used to generate e.g. a plot of pressure
% along certain wall.
%
% nodeIds = extractNodeIdsOnEdges(e, edgeIds)
%
% Input arguments:
% e - array of edge elements as defined in help to the importMeshGmsh
%     function.
% edgeIds - a row vector of edgeIds for which the node ids are to be
%     generated (e.g. ids of edges along which we later want to draw a
%     solution plot).
%
% Output arguments:
% nodeIds - global ids of mesh nodes lying on the edges specified by
%     edgeIds.
%
% Examples:
%       Draw pressure plot along edge 5, 6 and 7.
%           pressure = generatePressureData(u, p, t);
%           wallNodes = extractNodeIdsOnEdges(e, [5 6 7]);
%           scatter(p(1,wallNodes), pressure(wallNodes)')
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: EXTRACTDATAALONGLINE, GENERATEPRESSUREDATA, IMPORTMESHGMSH.

if(size(e,1)==7) % First order mesh
    nnodes = max(max(e(1:2,:)));
    nodeIds = zeros(1,nnodes);
    
    for faceId = edgeIds
        for edge = 1:size(e,2)
            if(e(5,edge)==faceId)
                dnodes = e([1 2],edge);
                
                nodeIds(dnodes) = 1;
            end
        end
    end
else % Second order mesh
    nnodes = max(max(e([1 2 8],:)));
    nodeIds = zeros(1,nnodes);
    
    for faceId = edgeIds
        for edge = 1:size(e,2)
            if(e(5,edge)==faceId)
                dnodes = e([1 2 8],edge);
                
                nodeIds(dnodes) = 1;
            end
        end
    end
end
    
nodeIds = find(nodeIds);

end