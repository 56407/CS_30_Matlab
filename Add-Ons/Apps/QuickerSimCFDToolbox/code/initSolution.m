function [u, convergence] = initSolution(p,t,velocityValue,pressureValue)

% initSolution - Set initial guess of pressure and velocity fields.
%
% This QuickerSim CFD Toolbox function is intended to init solution vector
% (both for velocity and pressure fields) with constant values for each
% field separately. The solver uses this initialization as a starting
% vector of the iterations. The function also initializes a convergence
% array which can be passed to the computeResiduals function.
%
% [u, convergence] = initSolution(p, t, velocityValue, pressureValue)
%
% Input arguments:
% p - array storing coordinates of all nodes in the mesh (detailed
%     description in help for importMeshGmsh function).
% t - array defining all elements of the mesh (detailed description also in
%     help for the importMeshGmsh function).
% velocityValue - a two (in 2-D case) or three (in 3-D case) element vector
%     with x-velocity value, y-velocity value and z-velocity value
%     (in 3-D case) which will be assigned as a starting approximation of
%     the velocity field in all nodes of the domain.
% pressureValue - a single value which will be assigned to all pressure
%     unknowns in the whole mesh.
%
% Output arguments:
% u - initialized solution vector of proper size with assigned values of
%     velocityValue and pressureValue. This vector can directly be passed
%     to the assembleNavierStokesMatrix as an initial approximation to the
%     convection velocity.
% convergence - an array which is intended to store during iteration
%     process ordinate number of the iteration and all residuals achieved
%     by the solver in each step. At the execution of initSolution this
%     array is set to [] and in this form can be directly passed as an
%     input argument to the computeResiduals function.
% 
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also ASSEMBLENAVIERSTOKESMATRIX2D, COMPUTERESIDUALS.

dim = size(p,1);
nnodes = size(p,2);

convergence = [];

if(dim == 2)
    % Define initial approximation
    nPnodes = max(max(t(1:3,:)));
    vx = velocityValue(1)*ones(nnodes,1);
    vy = velocityValue(2)*ones(nnodes,1);
    pressure = pressureValue*ones(nPnodes,1);
    u = [vx; vy; pressure];
else
    % Define initial approximation
    nPnodes = max(max(t(1:4,:)));
    vx = velocityValue(1)*ones(nnodes,1);
    vy = velocityValue(2)*ones(nnodes,1);
    vz = velocityValue(3)*ones(nnodes,1);
    pressure = pressureValue*ones(nPnodes,1);
    u = [vx; vy; vz; pressure];
end

end