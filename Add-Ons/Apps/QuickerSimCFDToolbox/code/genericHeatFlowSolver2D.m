% QuickerSim CFD Toolbox: Solve a stationary heat convection problem in 2-D

clc;
clear;

% Read and display mesh
[p,e,t] = importMeshGmsh('grids\heatExchanger.msh');
displayMesh2D(p,t);

% Convert mesh to second order
[p,e,t,nVnodes,nPnodes,indices] = convertMeshToSecondOrder(p,e,t);

% Define fluid properties
nu = 1e-6; % kinematic viscosity
lambda = 0.6; % thermal conductivity
rho = 1000; % fluid density
cp = 4190; % heat capacity
k = lambda/(rho*cp); % resulting thermal diffusivity

% Init solution and convergence criteria
[u, convergence] = initSolution(p,t,[0.05 0],0);
maxres = 1e-9;
maxiter = 25;

% Define inlet velocity profile
vel = [0.05 0];

% Iterate nonlinear terms
for iter = 1:maxiter
    % Assemble matrix and right hand side
    [NS, F] = assembleNavierStokesMatrix2D(p,e,t,nu,u(indices.indu),u(indices.indv),'nosupg');
    % [NS, F] = assembleNavierStokesMatrix2D(p,e,t,nu,u(indices.indu),u(indices.indv),'supgDoublyAsymptotic');
    
    % Apply boundary conditions
    [NS, F] = imposeCfdBoundaryCondition2D(p,e,t,NS,F,22,'inlet', vel);
    [NS, F] = imposeCfdBoundaryCondition2D(p,e,t,NS,F,[24 25],'slipAlongX');
    [NS, F] = imposeCfdBoundaryCondition2D(p,e,t,NS,F,[26 27 28 29 30],'wall');
    
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

% Display velocity field
figure(3);
displaySolution2D(p,t,u(indices.indu),'x-velocity');

% Solve heat transfer problem

% Assemble problem matrix
[D,F] = assembleDiffusionMatrix2D(p,t,k);
D = D + assembleScalarConvectionMatrix2D(p,t,k,u(indices.indu),u(indices.indv),'nosupg');

% Apply boundary conditions
Tinlet = 20;
Tpipe = 50;
[D,F] = imposeScalarBoundaryCondition2D(p,e,D,F,22,'value',Tinlet);
[D,F] = imposeScalarBoundaryCondition2D(p,e,D,F,[26 27 28 29 30],'value',Tpipe);

% Solve equations
T = D\F;

% Display solution
figure(4);
displaySolution2D(p,t,T,'Temperature field');