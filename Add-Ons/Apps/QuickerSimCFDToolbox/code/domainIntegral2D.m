function [J, av, area] = domainIntegral2D(p, t, u)

% domainIntegral2D - Compute integral of a scalar field over problem
%                    domain, average value and surface area of the whole
%                    domain.
%
% This QuickerSim CFD Toolbox function computes integral of a given scalar
% field u over whole domain and additionally returns average value of the
% field, as well as total surface area of the domain.
%
% [I, av, area] = domainIntegral2D(p, t, u);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% t - array of finite elements (see help of the importMeshGmsh function for
%     details).
% u - a column vector with scalar field to be integrated.
%
% Output arguments:
% I - a single scalar value of the integral over whole domain.
% av - average value of the u field on whole domain area.
% area - surface area of the problem domain.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: BOUNDARYFLUX2D, BOUNDARYINTEGRAL2D, CONVERTMESHTOSECONDORDER,
%           IMPORTMESHGMSH.


nelements = size(t,2);
S = zeros(1,nelements);
mv = zeros(1,nelements);

for el = 1:nelements
    node = t(1:3,el);
    v1 = [p(:,node(2))-p(:,node(1)); 0];
    v2 = [p(:,node(3))-p(:,node(1)); 0];
    
    S(el) = norm(cross(v1,v2));
    mv(el) = sum(u(node))/3;
end

S = 0.5*S;

J = sum(mv.*S);
area = sum(S);
av = J/area;

end