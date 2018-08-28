function [Ftotal, Fviscous] = computeForce(p, e, t, u, nu, wallIds)

% computeForce - Compute force acting from fluid on given surfaces.
%
% This QuickerSim CFD Toolbox function calculates both total and viscous
% force acting from fluid on the boundary given by wallIds. The results
% should be multiplied with fluid density to obtain force in Newtons.
%
% [Ftotal, Fviscous] = computeForce(p, e, t, u, nu, wallIds);
%
% Input arguments:
% p, e, t - arrays which store computational mesh (all details of  the PET
%     mesh format described in help for importMeshGmsh function).
% u - a column vector with solution to the velocity and pressure field.
% nu - kinematic viscosity value in one of two formats:
%     1. a scalar value constant in the whole domain,
%     2. a vector of nu values for each node in the mesh.
% wallIds - row vector of ids of edges forming the boundary on which the
%     force is to be computed.
% 
% Output arguments:
% Ftotal - a row vector of the resulting force with x, y and z component
%          (in 3-D).
% Fviscous - a row vector with x, y and z (if in 3-D) component of the
%            viscous force term. Hence, the net pressure part of force can
%            be computed as Ftotal-Fviscous.
%
% Examples:
%       Compute force acting on wall composed from boundaries 5 and 7.
%       F = computeForce(p, e, t, u, nu, [5 7]);
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: BOUNDARYINTEGRAL2D, BOUNDARYFLUX2D, COMPUTEMOMENT,
%           IMPORTMESHGMSH.

dim = size(p,1);

if(dim == 2)
    % Compute pressure force
    pressure = generatePressureData(u, p, t);
    Fpressure = boundaryIntegral2D(p,e,pressure,wallIds,'normalVector');
    %disp(Fpressure);

    % Compute viscous force
    nnodes = size(p,2);
    ugrad = solutionGradient2D(p,t,u(1:nnodes));
    vgrad = solutionGradient2D(p,t,u(nnodes+1:2*nnodes));
    
    % -p*n+2*nu*D*n
    D11 = ugrad(:,1);
    D12 = 0.5*(ugrad(:,2)+vgrad(:,1));
    D21 = D12;
    D22 = vgrad(:,2);

    if(length(nu)==1)
        nu = repmat(nu,nnodes,1);
    end

%     f1 = [-nu.*ugrad(:,1); -nu.*ugrad(:,2)];
%     f2 = [-nu.*vgrad(:,1); -nu.*vgrad(:,2)];

    f1 = [-2*nu.*D11; -2*nu.*D12];
    f2 = [-2*nu.*D21; -2*nu.*D22];

    fx_viscous = boundaryFlux2D(p,e,f1,wallIds);
    fy_viscous = boundaryFlux2D(p,e,f2,wallIds);

    Fviscous = [fx_viscous fy_viscous];
    Ftotal = Fpressure + Fviscous;
else % 3-D
    % Compute pressure force
    pressure = generatePressureData(u,p,t);
    Fpressure = boundaryIntegral(p,e,pressure,wallIds,'normalVector');
    
    % Compute viscous force from n*(2*mu*D)
    nnodes = size(p,2);
    ugrad = solutionGradient(p,t,u(1:nnodes));
    vgrad = solutionGradient(p,t,u((nnodes+1):(2*nnodes)));
    wgrad = solutionGradient(p,t,u((2*nnodes+1):(3*nnodes)));
    
    % -p*n + n*2*mu*D
    D11 = ugrad(:,1);
    D12 = 0.5*(ugrad(:,2)+vgrad(:,1));
    D13 = 0.5*(ugrad(:,3)+wgrad(:,1));
    D21 = D12;
    D22 = vgrad(:,2);
    D23 = 0.5*(vgrad(:,3)+wgrad(:,2));
    D31 = D13;
    D32 = D23;
    D33 = wgrad(:,3);
    
    if(length(nu)==1)
        nu = repmat(nu,nnodes,1);
    end
    
    f1 = -[2*nu.*D11; 2*nu.*D12; 2*nu.*D13];
    f2 = -[2*nu.*D21; 2*nu.*D22; 2*nu.*D23];
    f3 = -[2*nu.*D31; 2*nu.*D32; 2*nu.*D33];
    
    fx_v = boundaryFlux(p,e,f1,wallIds);
    fy_v = boundaryFlux(p,e,f2,wallIds);
    fz_v = boundaryFlux(p,e,f3,wallIds);
    
    Fviscous = [fx_v fy_v fz_v];
    Ftotal = Fpressure+Fviscous;
end

end