function [p, e, t, nVnodes, nPnodes, indices] = convertMeshToSecondOrder(p, e, t)

% convertMeshToSecondOrder - Convert linear grid to second order grid.
%
% This QuickerSim CFD Toolbox function converts a first order linear mesh
% read from an external file (e.g. Gmsh etc.) to the second order mesh
% which is needed for fluid flow simulation with P2-P1 elements (so called
% Taylor-Hood pair) with second order shape functions for velocity
% discretization and linear shape functions for pressure approximation.
%
% [p, e, t, nVnodes, nPnodes, indices] = convertMeshToSecondOrder(p, e, t)
%
% Input arguments:
% p - array of nodal coordinates in the linear mesh.
% e - array of boundary elements in the linear mesh.
% t - array of linear domain elements (triangles/tetrahedra).
%
% Output arguments:
% p, e, t - arrays defining mesh - see help for importMeshGmsh function for
%     detailed description of mesh PET format.
% p - as above with additional second order nodes in the middle of each
%     edge of every element.
% e - as above with added indices of new nodes. For a linear 2-D mesh this
%     array has 7 rows, for a second order grid the number of rows equals
%     8 - the last row stores the indices of the new nodes added in the
%     middle of each edge. For the 3-D case number of rows in e is either 6
%     or 9 for linear and second order grid, respectively. For details of
%     node numbering refer to documentation.
% t - as above. For a 2-D linear mesh this array stores indices of element
%     nodes in rows 1 to 3 and the id of the domain in the last 4-th row. In
%     case of second order grid the corner nodes of each triangle are still
%     stored in rows 1 to 3 and rows 4 to 6 store correspondingly indices of
%     new second order nodes in the following order: node 4 (placed between
%     original 1st and 2nd node), node 5 (placed between original 2nd and 3rd
%     node) and node 6 (placed between original 3rd and 1st node). In 3-D
%     number of rows in t equals 5 for a linear mesh (4 element nodes and
%     domain id) or 11 for second order mesh (10 ids of element nodes and
%     domain id in the last row). For details of node numbering refer to
%     documentation.
% nVnodes - total number of nodes associated with velocity unknowns (so called
%     velocity nodes - i.e. all nodes in the second order mesh).
% nPnodes - total number of nodes associated with pressure unknowns (so
%     called pressure nodes - i.e. all nodes of the original linear grid).
% indices - a structure containing for each node indices of the x-velocity,
%     y-velocity, z-velocity (only for 3-D case) and pressure unknown 
%     in the global solution vector u.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also IMPORTMESHGMSH, GENERATEINDICES2D.

dim = size(p,1);

nnodes = size(p,2);
nelements = size(t,2);

Edges = cell(nnodes,1);
AddNodes = cell(nnodes,1);

for i = 1:nnodes
    Edges{i} = i;
    AddNodes{i} = i;
end

if(dim == 2)
    t = [t; zeros(3, nelements)];
    t(7,:) = t(4,:);
else % dim == 3
    t = [t; zeros(6, nelements)];
    t(11,:) = t(5,:);
end
el = nnodes+1;

if(dim == 2)
    for i = 1:nelements
        for j = 1:3
            n1 = t(j,i);
            n2 = t(mod(j,3)+1,i);
            nloc = sort([n1 n2]);
            n1 = nloc(1);
            n2 = nloc(2);

            ind = find(Edges{n1}==n2);

            if(ind)
                t(3+j,i) = AddNodes{n1}(ind);
            else
                Edges{n1}(end+1) = n2;
                AddNodes{n1}(end+1) = el;
                t(3+j,i) = el;
                el = el+1;
            end
        end
    end
else % dim == 3
    edgeNodes = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4];
    for i = 1:nelements
        for j = 1:6
            n1 = t(edgeNodes(j,1),i);
            n2 = t(edgeNodes(j,2),i);
            nloc = sort([n1 n2]);
            n1 = nloc(1);
            n2 = nloc(2);

            ind = find(Edges{n1}==n2);

            if(ind)
                t(4+j,i) = AddNodes{n1}(ind);
            else
                Edges{n1}(end+1) = n2;
                AddNodes{n1}(end+1) = el;
                t(4+j,i) = el;
                el = el+1;
            end
        end
    end
end

el = el-1;

przyrost = el-nnodes;
nnodesfull = el;

if(dim == 2)
    p = [p zeros(2,przyrost)];
    
    for i = 1:nelements
        for j = 1:3
            n1 = t(j,i);
            n2 = t(mod(j,3)+1,i);

            p(:,t(j+3,i)) = 0.5*(p(:,n1)+p(:,n2));
        end
    end
else
    p = [p zeros(3,przyrost)];
    
    edgeNodes = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4];
    for i = 1:nelements
        for j = 1:6
            n1 = t(edgeNodes(j,1),i);
            n2 = t(edgeNodes(j,2),i);

            p(:,t(j+4,i)) = 0.5*(p(:,n1)+p(:,n2));
        end
    end
end


% Enrich edges
if(dim == 2)
    nEdges = size(e,2);
    e = [e; zeros(1,nEdges)];

    for i = 1:nEdges
        n1 = e(1,i);
        n2 = e(2,i);
        nloc = sort([n1 n2]);
        n1 = nloc(1);
        n2 = nloc(2);
        ind = find(Edges{n1}==n2);
        e(8,i) = AddNodes{n1}(ind);
    end

    nVnodes = size(p,2);
    nPnodes = max(max(t(1:3,:)));
    indices = generateIndices2D(p,t);
else % dim == 3
    nEdges = size(e,2);
    e = [e(1:3,:); zeros(3,nEdges); e(4:6,:)];

    for i = 1:nEdges
        % Popraw cialo tej petli
        n1 = e(1,i);
        n2 = e(2,i);
        nloc = sort([n1 n2]);
        n1 = nloc(1);
        n2 = nloc(2);
        e(4,i) = AddNodes{n1}(Edges{n1}==n2);
        
        n1 = e(2,i);
        n2 = e(3,i);
        nloc = sort([n1 n2]);
        n1 = nloc(1);
        n2 = nloc(2);
        e(5,i) = AddNodes{n1}(Edges{n1}==n2);
        
        n1 = e(1,i);
        n2 = e(3,i);
        nloc = sort([n1 n2]);
        n1 = nloc(1);
        n2 = nloc(2);
        e(6,i) = AddNodes{n1}(Edges{n1}==n2);
    end

    nVnodes = size(p,2);
    nPnodes = max(max(t(1:4,:)));
    indices = generateIndices(p,t);
end

end