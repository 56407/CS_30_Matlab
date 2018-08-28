% QuickerSim CFD Toolbox: Solve a stationary Navier-Stokes problem in 2-D

clc;
clear;

% Read and display mesh
[p,e,t] = importMeshGmsh('grids/naca5418.msh');
displayMesh2D(p,t);

% Convert mesh to second order
[p,e,t,nVnodes,nPnodes,indices] = convertMeshToSecondOrder(p,e,t);

% Define fluid
nu = 0.02;

% Init solution and convergence criteria
[u, convergence] = initSolution(p,t,[1 0],0);
maxres = 1e-3;
maxiter = 25;

% Define inlet velocity profile
vel = [1 0];

% Iterate nonlinear terms
for iter = 1:maxiter
    % Compute shear rate and viscosity for non-Newtonian fluids if needed
    % gamma = shearRate(p,t,u);
    % nu = nuinf + (nu0-nuinf)*(1+(lambda*gamma).^2).^((n-1)/2);
    
    % Assemble matrix and right hand side
    [NS, F] = assembleNavierStokesMatrix2D(p,e,t,nu,u(indices.indu),u(indices.indv),'nosupg');
    % [NS, F] = assembleNavierStokesMatrix2D(p,e,t,nu,u(indices.indu),u(indices.indv),'supgDoublyAsymptotic');
    
    % Apply boundary conditions
    [NS, F] = imposeCfdBoundaryCondition2D(p,e,t,NS,F,10,'inlet', [1 0]);
    [NS, F] = imposeCfdBoundaryCondition2D(p,e,t,NS,F,12,'slipAlongX');
    [NS, F] = imposeCfdBoundaryCondition2D(p,e,t,NS,F,13,'wall');
    
    % Compute and plot residuals
    [stop, convergence] = computeResiduals(NS,F,u,size(p),convergence,maxres);
    plotResiduals(convergence,2);
    
    % Break if solution converged
    if(stop)
        break;
    end
    
    % Solve equations
    u = NS\F;
end

% Plot solution
figure(3);
displaySolution2D(p,t,u(indices.indu),'x-velocity');