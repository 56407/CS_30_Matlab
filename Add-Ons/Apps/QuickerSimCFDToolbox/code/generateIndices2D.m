function [indices] = generateIndices2D(p, t)

% generateIndices2D - Generate global ids of pressure and velocity unknowns
%
% This QuickerSim CFD Toolbox function generates indices for easy
% manipulation of the solution vector u.
%
% indices = generateIndices2D(p,t);
%
% Input arguments:
% p - array of nodal coordinates (for detailed description see help of the
%     importMeshGmsh function).
% t - array of finite elements (for detailed description see help for the
%     importMeshGmsh function).
% 
% Output arguments:
% indices - a Matlab structure with fiels:
%           indices.indu - with indices of x-velocity unknowns,
%           indices.indv - with indices of y-velocity unknowns,
%           indices.indp - with indices of pressure unknowns.
% The general structure of the solution vector u in this Toolbox is that
% the first 1 to nVnodes (where nVnodes stands for total number of velocity
% nodes - see help for convertMeshToSecondOrder for detailed description of
% nVnodes) elements in the solution vector u correspond to x-velocity
% values in the subsequent nodes of the mesh, the next elements in u, this
% means u(nVnodes+1:2*nVnodes) contain values of the y-velocity for each
% mesh node and then the last elements, this means u(2*nVnodes+1:end)
% contain pressure values for each mesh node (note here that only first
% order nodes contain pressure solution - do not refer to mid-edge nodes
% when requesting pressure values, since this will cause an error with
% false indices in indices.indp vector). Function generatePressureData
% might be helpful.
%
% Examples:
% 1. Extract y-velocity value for the 179th node in the mesh:
%       indices = generateIndices2D(p,t);
%       v179 = u(indices.indv(179));
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: CONVERTMESHTOSECONDORDER, EXPORTTOGMSH2D, GENERATEPRESSUREDATA,
%           IMPORTMESHGMSH.

% Pressure nodes
pNodes = t(1:3,:);

% pNodesSize = size(pNodes);
% pNodesList = sort(reshape(pNodes,1,pNodesSize(1)*pNodesSize(2)));
% pNodesUniqueList = unique(pNodesList);

vNodesNumber = size(p,2);
pNodesNumber = max(max(pNodes));

% disp([pNodesNumber, vNodesNumber])

indu = 1:vNodesNumber;
indv = (vNodesNumber+1):(2*vNodesNumber);
%indw = (2*vNodesNumber+1):(3*vNodesNumber);
%indp = zeros(1,vNodesNumber);
indp = zeros(1,pNodesNumber);
indp(1:pNodesNumber) = (1:pNodesNumber)+2*vNodesNumber;

indices = struct('indp',indp,'indu',indu,'indv',indv);

end