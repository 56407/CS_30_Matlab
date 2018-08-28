% QuickerSim CFD Toolbox: Solve a stationary heat conduction problem in 2-D

clc;
clear;

% Import mesh
[p,e,t] = importMeshGmsh('grids/rectWithHoles.msh');

% Define convergence criteria
maxiter = 10;
maxres = 1e-3;

% Init solution with T = 20
[T,convergence] = initScalarSolution(p,20);

% Start iteration loop to iterate nonlinearity
for iter = 1:maxiter
    % Compute thermal conductivity
    lambda = 20+0.4*T; % W/(mK)

    % Assemble problem matrix
    [D,F] = assembleDiffusionMatrix2D(p,t,lambda);

    % Apply boundary conditions
    [D,F] = imposeScalarBoundaryCondition2D(p,e,D,F,18,'value',20);
    [D,F] = imposeScalarBoundaryCondition2D(p,e,D,F,17,'flux',1000);
    
    % Compute and plot residuals
    [stop, convergence] = computeScalarResiduals(D,F,T,convergence,maxres);
    plotScalarResiduals(convergence,2);
    
    % Break if solution converged
    if(stop)
        break;
    end

    % Solve linearized equations
    T = D\F;
end

figure(4);
displaySolution2D(p,t,T,'Temperature [deg C]');