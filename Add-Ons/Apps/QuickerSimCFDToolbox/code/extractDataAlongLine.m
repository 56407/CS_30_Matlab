function [uValues, linep, linee] = extractDataAlongLine(p, t, u, lineX, lineY, npoints, edgeId)

% extractDataAlongLine - Extract field values at arbitrary point in the
%                        mesh.
%
% This QuickerSim CFD Toolbox function is intended mainly to extract data
% along any arbitrary line for convenient plot creation of flow profile at
% any given location.
%
% [v,lp,le] = extractDataAlongLine(p,t,u,lineX,lineY,npoints,edgeId);
%
% Input arguments:
% p - array of nodal coordinates described in details in help to the
%     importMeshGmsh function.
% t - array of triangular finite elements described exactly in help to the
%     importMeshGmsh function.
% u - a nnodes-element column vector representing any scalar field (nnodes
%     stands for the total number of mesh nodes). Each element of the u
%     vector represents value of the u scalar field at this particular mesh
%     node.
% lineX - a 2-element row vector with x coordinates of the beginning and
%     end of the line along which we want to distribute points and extract
%     data.
% lineY - a 2-element row vector with y coordinates of the beginning and
%     end of the line along which we want to distribute points and extract
%     data.
% npoints - number of points which shall be distributed along line and in
%     which we want to get u-field values.
% edgeId - an arbitrary number (an id) which we want to be assigned to the
%     newly created edge domain. This issue is explained more clearly in
%     the 'Output arguments' section.
%
% Output arguments:
% v - a row vector of extracted values of u-field along the line.
% lp - a 2-by-npoints array of coordinates of locations at which the
%     u-field data have been extracted. The format of this array is
%     identical with the format of p (described in help to importMeshGmsh
%     function) except for the fact that it does not store mesh nodes, but
%     the points along extraction line.
% le - a 7-by-(npoints-1) array of edge definition. Each line segment
%     between two extraction points is defined as an edge element and
%     placed in le array (format is identical with the format of e array
%     described in help to the importMeshGmsh function). Rows 1 and 2 are
%     filled in with ids of newly created extraction points (whose
%     coordinates are accessible through lp array), rows 3 and 4 contain
%     nonvalid, arbitrary data, row 5 contains edgeId value, rows 6 and 7
%     are assigned with 1 (subdomain label which is always assumed one in
%     the Lite version of the Toolbox). Defining lp and le arrays enables
%     computing boundary integrals with the boundaryFlux2D or
%     boundaryIntegral2D functions, thus allowing for more advanced
%     analysis of results. For computing flux values and other orientation
%     dependent quantities note that outward unit normal vector on the edge
%     is always defined as pointing to the left when marching along the
%     edge (line) from its beginning to its end.
% 
% Examples:
%       For valuable practical examples of usage of this function refer to
%       courses and tutorials.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: BOUNDARYFLUX2D, BOUNDARYINTEGRAL2D, EXTRACTNODEIDSONEDGES,
%           IMPORTMESHGMSH.

x = linspace(lineX(1), lineX(2), npoints);
y = linspace(lineY(1), lineY(2), npoints);
uValues = zeros(1,length(x));

for i = 1:length(x)
    % Znajdz trojkat
    [el, b1, b2, b3] = pointInTriangle(p,t,[x(i) y(i)]);

    if(el==0)
        uValues(i) = nan;
    else
        uValues(i) = b1*u(t(1,el)) + b2*u(t(2,el)) + b3*u(t(3,el));
        %uValues(i) = (u(t(1,el))+u(t(2,el))+u(t(3,el)))/3;
    end
end

% Stworz tablice e i p
linee = zeros(7,npoints-1);

linee(1,:) = 1:npoints-1;
linee(2,:) = 2:npoints;
linee(5,:) = edgeId;
linee(6:7,:) = 1;

linep = [x; y];


end