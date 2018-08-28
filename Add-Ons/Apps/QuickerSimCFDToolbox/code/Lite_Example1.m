function [] = Lite_Example1()

clc;
clear;

% Read mesh
load backwardFacingStepGeometry2D.mat;

% Convert mesh to second order
[p,e,t,nVnodes,nPnodes,indices] = convertMeshToSecondOrder(p,e,t);

% Define fluid
nu = 0.02;

% Initialize solution and set convergence criteria
[u, convergence] = initSolution(p,t,[0.5 0],0);
maxres = 1e-3;
maxiter = 25;

% Define inlet velocity profile
vel = @(x,y)([6*(y)*(1-y) 0]);

% Iterate nonlinear terms
for iter = 1:maxiter
    % Assemble matrix and right hand side
    [NS, F] = assembleNavierStokesMatrix2D(p,e,t,nu,u(indices.indu),u(indices.indv),'nosupg');
    
    % Apply boundary conditions
    [NS, F] = imposeCfdBoundaryCondition2D(p,e,t,NS,F,[1 3 5 6 7],'wall');
    [NS, F] = imposeCfdBoundaryCondition2D(p,e,t,NS,F,2,'inlet', vel);
    
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
figure(1)
displaySolution3D(p,t,u(indices.indu),'x-velocity');
xlim([-1 7])
ylim([-3 3])

% Export results
exportToGmsh2D('example1_result.msh',u(indices.indu),p,t,'x-velocity');


disp('Congratulations! You can view the source code by typing:');
disp('edit Lite_Example1');
disp('or go directly to www.quickersim.com/cfd-toolbox-for-matlab/index');
disp('to download your first geometry creation, meshing, solving and');
disp('postprocessing tutorial and some additonal source code examples');
disp('that will teach you using this QuickerSim CFD Toolbox.');