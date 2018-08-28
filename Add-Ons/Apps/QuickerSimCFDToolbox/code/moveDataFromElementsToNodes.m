function f = moveDataFromElementsToNodes(p,t,u)

% moveDataFromElementsToNodes - Interpolate data values stored in element
%                               centroids to nodal locations.
%
% This QuickerSim CFD Toolbox function enables interpolating data from
% element cetroids and store them at nodal locations for a 2-D mesh.
%
% f = moveDataFromElementsToNodes(p,t,u);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% t - array of finite elements (see help of the importMeshGmsh function for
%     details).
% u - vector containing field values expressed at element centroids.
%
% Output arguments:
% f - vector containing field values expressed at nodes
% 
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: IMPORTMESHGMSH, IMPOSESCALARBOUNDARYCONDITION2D,
%           SOLUTIONGRADIENT2D.

nnodes = size(p,2);
f = zeros(nnodes,size(u,2));
s = zeros(nnodes,1);

% First order mesh
if(size(t,1)==4)
    for el = 1:size(t,2)
        elementNodes = t(1:3,el);
        f(elementNodes,:) = f(elementNodes,:) + u(el,:);
        s(elementNodes) = s(elementNodes)+1;
    end
else
    for el = 1:size(t,2)
        elementNodes = t(1:6,el);
        f(elementNodes,:) = f(elementNodes,:) + u(el,:);
        s(elementNodes) = s(elementNodes)+1;
    end
end

f = f./repmat(s,1,size(u,2));

end