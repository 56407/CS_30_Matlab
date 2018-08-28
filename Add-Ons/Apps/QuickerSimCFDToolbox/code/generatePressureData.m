function [pressure] = generatePressureData(u, p, t)

% generatePressureData - Compute pressure in all mesh nodes.
%
% This QuickerSim CFD Toolbox function interpolates data from pressure
% nodes to all nodes in the mesh. Simulations in this toolbox are performed
% with a Hood-Taylor pair of elements, fulfilling the LBB (inf-sup)
% condition with second order shape functions for velocity and linear
% shape functions for pressure. This means that we only get computed
% pressure values in the three corner nodes of a triangle or four corner
% nodes of tetrahedron. This function interpolates these values to the
% second order nodes in the middle of each edge by means of linear
% interpolation. The output of this function can directly be used in the
% call to the exportToGmsh2D or exportToGmsh function.
%
% pressure = generatePressureData(u, p, t);
%
% Input arguments:
% u - a dim*nVnodes+nPnodes-by-1 solution vector in the usual CFD Toolbox
%     form, where dim is 2 for 2-D case and 3 for 3-D case (for a more
%     detailed discussion of nVnodes and nPnodes see help for
%     convertMeshToSecondOrder function).
% p - array of nodal coordinates (details of the format in help for
%     importMeshGmsh function).
% t - array of elements (details of the format in help for importMeshGmsh
%     function).
% 
% Output arguments:
% pressure - a nVnodes-by-1 vector with pressure values calculated for all
%     nodes in the mesh.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: CONVERTMESHTOSECONDORDER, EXPORTTOGMSH2D, IMPORTMESHGMSH.

nVnodes = size(p,2);

pressure = zeros(nVnodes,1);

nelements = size(t,2);

dim = size(p,1);

if(dim == 3)
    for el = 1:nelements
        elementNodes = t(1:10,el);

        pressure(elementNodes(1:4)) = u(3*nVnodes+elementNodes(1:4));
        pressure(elementNodes(5)) = 0.5*(pressure(elementNodes(1))+pressure(elementNodes(2)));
        pressure(elementNodes(6)) = 0.5*(pressure(elementNodes(2))+pressure(elementNodes(3)));
        pressure(elementNodes(7)) = 0.5*(pressure(elementNodes(1))+pressure(elementNodes(3)));
        pressure(elementNodes(8)) = 0.5*(pressure(elementNodes(1))+pressure(elementNodes(4)));
        pressure(elementNodes(9)) = 0.5*(pressure(elementNodes(2))+pressure(elementNodes(4)));
        pressure(elementNodes(10)) = 0.5*(pressure(elementNodes(3))+pressure(elementNodes(4)));
    end
elseif(dim == 2)
    for el = 1:nelements
        elementNodes = t(1:6,el);

        pressure(elementNodes(1:3)) = u(2*nVnodes+elementNodes(1:3));
        pressure(elementNodes(4)) = 0.5*(pressure(elementNodes(1))+pressure(elementNodes(2)));
        pressure(elementNodes(5)) = 0.5*(pressure(elementNodes(2))+pressure(elementNodes(3)));
        pressure(elementNodes(6)) = 0.5*(pressure(elementNodes(1))+pressure(elementNodes(3)));
    end
end

end