function [omega, omegamag] = vorticity(p,t,u)

% vorticity - Compute vorticity.
%
% This QuickerSim CFD Toolbox function computes vorticity vector and
% magnitude for a given flow field.
%
% [omega, omegamag] = vorticity(p, t, u)
%
% Input arguments:
% p - array of nodal coordinates described in details in help to the
%     importMeshGmsh function.
% t - array of triangular finite elements described exactly in help to the
%     importMeshGmsh function.
% u - solution vector to the velocity field.
%
% Output arguments:
% omega - nVnodes-by-1 element vector with z-vorticity values for each node
%         (in 2-D case) or nVnodes-by-3 element array with x, y and z
%         component of vorticity for each node (in case of a 3-D case).
%         omega is calculated as rot(u).
% omegamag - magnitude of vorticity calculated for each node from:
%            omegamag = sqrt(2*Omega_ij*Omega_ij), where Omega_ij is given
%            by: Omega_ij = 1/2*(d_ui/d_xj - d_uj/d_xi)
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: COMPUTEFORCE, COMPUTEMOMENT, EXPORTTOGMSH2D, IMPORTMESHGMSH,
%           SOLUTIONGRADIENT2D, SHEARRATE.

dim = size(p,1);
nnodes = size(p,2);

if(dim == 2)
    ugrad = solutionGradient2D(p,t,u(1:nnodes));
    vgrad = solutionGradient2D(p,t,u((nnodes+1):(2*nnodes)));

    %ux = ugrad(:,1);
    uy = ugrad(:,2);
    vx = vgrad(:,1);
    %vy = vgrad(:,2);

    %omega11 = 0.5*(ux-ux);
    omega21 = 0.5*(vx-uy);
    omega12 = 0.5*(uy-vx);
    %omega22 = 0.5*(vy-vy);

    omegamag = sqrt(2*omega21.^2 + 2*omega12.^2);
    omega = vx-uy;
else % dim == 3
    ugrad = solutionGradient(p,t,u(1:nnodes));
    vgrad = solutionGradient(p,t,u((nnodes+1):(2*nnodes)));
    wgrad = solutionGradient(p,t,u(2*(nnodes+1):(3*nnodes)));
    
    wy = wgrad(:,2);
    vz = vgrad(:,3);
    uz = ugrad(:,3);
    wx = wgrad(:,1);
    vx = vgrad(:,1);
    uy = ugrad(:,2);
    
    omega12 = 0.5*(uy-vx);
    omega13 = 0.5*(uz-wx);
    omega21 = -omega12;
    omega23 = 0.5*(vz-wy);
    omega31 = -omega13;
    omega32 = -omega23;
    
    omegax = wy-vz;
    omegay = uz-wx;
    omegaz = vx-uy;
    
    omega = [omegax omegay omegaz];
    omegamag = sqrt(2*(omega12.^2 + omega13.^2 + omega21.^2 + omega23.^2 + omega31.^2 + omega32.^2));
end

end