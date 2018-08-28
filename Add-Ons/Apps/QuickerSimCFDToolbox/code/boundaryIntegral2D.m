function [J, av] = boundaryIntegral2D(p,e,u,edgeIds,varargin)

% boundaryIntegral2D - Compute integral and average value of a scalar field
%                      at boundary.
%
% This QuickerSim CFD Toolbox function computes boundary integral of a
% scalar field over specified edges. In addition it also returns average
% value of that quantity over the specified boundary. It can also compute
% integral of the u field times the boundary unit normal vector.
%
% [I, av] = boundaryIntegral2D(p, e, u, edgeIds);
% [I, av] = boundaryIntegral2D(p, e, u, edgeIds, integralType);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% e - array of edge elements lying on the boundary of the mesh (see help of
%     the importMeshGmsh function for details).
% u - a column vector with scalar field to be integrated.
% edgeIds - a row vector of ids of the edges over which integration shall
%     be performed.
% integralType - a string that specifies type of the boundary integral to
%     be computed. Currently, the only option which is supported is
%     'normalVector'. In that case, the function computes boundary integral
%     of the scalar field u times the outward unit normal vector and
%     returns a two-element vector result.
%
% Output arguments:
% I - a single scalar value of the boundary integral over given edges.
%     If integralType == 'normalVector' I is a 1-by-2 vector with values of
%     boundary integral of u times nx and u times ny on corresponding
%     components.
% av - average value of the u field on given boundary part.
%     If integralType == 'normalVector' av is a 1-by-2 vector of I divided
%     by the total boundary length.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: BOUNDARYFLUX2D, CONVERTMESHTOSECONDORDER, DOMAININTEGRAL2D,
%           IMPORTMESHGMSH.

if(~isempty(varargin))
    % disp('Licz z normalna!');
    nedges = size(e,2);
    J = [0 0];
    i = 0;
    Ltotal = 0;

    for edgeId = edgeIds
        for edge = 1:nedges
            if(e(5,edge) == edgeId)
                node = e(1:2,edge);
                v = p(:,node(2))-p(:,node(1));
                nx = v(2);
                ny = -v(1);
                L = norm(v);
                Ltotal = Ltotal + L;
                nx = nx/L;
                ny = ny/L;
                %J = J + 0.5*sum(u(node))*L;
                utmp = 0.5*sum(u(node))*L;
                J = J + [utmp*nx, utmp*ny];
                i = i+1;
            end
        end
    end
    av = J/Ltotal;
else
    nedges = size(e,2);
    J = 0;
    i = 0;
    Ltotal = 0;

    for edgeId = edgeIds
        for edge = 1:nedges
            if(e(5,edge) == edgeId)
                node = e(1:2,edge);
                v = p(:,node(2))-p(:,node(1));
                L = norm(v);
                Ltotal = Ltotal+L;
                J = J + 0.5*sum(u(node))*L;
                i = i+1;
            end
        end
    end
    av = J/Ltotal;
end

end