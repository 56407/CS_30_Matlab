function gamma = shearRate(p,t,u)

% shearRate - Compute shear rate.
%
% This QuickerSim CFD Toolbox function computes shear rate from given flow
% field.
%
% gamma = shearRate(p, t, u)
%
% Input arguments:
% p - array of nodal coordinates described in details in help to the
%     importMeshGmsh function.
% t - array of finite elements.
% u - solution vector to the velocity field.
%
% Output arguments:
% gamma - a nnodes-by-1 element vector with values of shear rate at every
%     node in the mesh. Shear rate is defined as:
%           gamma = sqrt(2*gamma_ij*gamma_ij), where:
%           gamma_ij = 0.5*(d/dx_j(u_i) + d/dx_i(u_j))
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: COMPUTEFORCE, COMPUTEMOMENT, EXPORTTOGMSH2D, IMPORTMESHGMSH,
%           SOLUTIONGRADIENT2D, VORTICITY.

dim = size(p,1);
nnodes = size(p,2);

if(dim == 2)
    ugrad = solutionGradient2D(p,t,u(1:nnodes));
    vgrad = solutionGradient2D(p,t,u((nnodes+1):(2*nnodes)));

    ux = ugrad(:,1);
    uy = ugrad(:,2);
    vx = vgrad(:,1);
    vy = vgrad(:,2);

    gamma11 = ux;
    gamma12 = 0.5*(uy+vx);
    gamma21 = gamma12;
    gamma22 = vy;

    gamma = sqrt(2*gamma11.^2 + 2*gamma12.^2 + 2*gamma21.^2 + 2*gamma22.^2);
else % dim = 3
    ugrad = solutionGradient(p,t,u(1:nnodes));
    vgrad = solutionGradient(p,t,u((nnodes+1):(2*nnodes)));
    wgrad = solutionGradient(p,t,u((2*nnodes+1):(3*nnodes)));
    
    ux = ugrad(:,1);
    uy = ugrad(:,2);
    uz = ugrad(:,3);
    vx = vgrad(:,1);
    vy = vgrad(:,2);
    vz = vgrad(:,3);
    wx = wgrad(:,1);
    wy = wgrad(:,2);
    wz = wgrad(:,3);
    
    e11 = ux;
    e12 = 0.5*(uy+vx);
    e13 = 0.5*(uz+wx);
    e21 = e12;
    e22 = vy;
    e23 = 0.5*(vz+wy);
    e31 = e13;
    e32 = e23;
    e33 = wz;
    
    gamma = sqrt(2*(e11.^2+e12.^2+e13.^2+e21.^2+e22.^2+e23.^2+e31.^2+e32.^2+e33.^2));
end

end