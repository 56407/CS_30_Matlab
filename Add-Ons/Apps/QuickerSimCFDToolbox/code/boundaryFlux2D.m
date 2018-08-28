function [J] = boundaryFlux2D(p,e,u,edgeIds)

% boundaryFlux2D - Compute flux of a vector field through boundary.
%
% This QuickerSim CFD Toolbox is used in postprocessing to compute the net
% flux of a vector quantity through domain boundary. This function computes
% a boundary integral of the term u_x * nx + u_y * ny, where nx and ny are
% the outward unit normal vector components and u_x and u_y the components
% of the u vector field.
%
% I = boundaryFlux2D(p,e,u,egdeIds);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% e - array of edge elements lying on the boundary of the mesh (see help of
%     the importMeshGmsh function for details).
% u - a column vector with 2*nnodes, where nnodes denotes total number of
%     nodes in the mesh. This vector includes values of the vector field to
%     be integrated. The first nnodes elements of the u vector must include
%     x-components of the field and the following nnodes the y-component
%     values.
% edgeIds - a row vector of ids of the edges over which integration shall
%     be performed.
%
% Output arguments:
% I - a single scalar value of the flux through the given boundary.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: BOUNDARYINTEGRAL2D, CONVERTMESHTOSECONDORDER, DOMAININTEGRAL2D,
%           IMPORTMESHGMSH.

nedges = size(e,2);
nnodes = size(p,2);
J = 0;
i = 0;

for edgeId = edgeIds
    for edge = 1:nedges
        if(e(5,edge) == edgeId)
            node = e(1:2,edge);
            v = p(:,node(2))-p(:,node(1));
            nx = v(2);
            ny = -v(1);
            L = norm(v);
            nx = nx/L;
            ny = ny/L;
            %J = J + 0.5*sum(u(node))*L;
            ux = 0.5*sum(u(node))*L;
            uy = 0.5*sum(u(nnodes+node))*L;
            J = J + ux*nx + uy*ny;
            i = i+1;
        end
    end
end

end