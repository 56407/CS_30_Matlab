% QuickerSim CFD Toolbox: Solve flow in porous media with Darcy's law

% Import mesh and convert to second order mesh
[p,e,t] = importMeshGmsh('grids/porous.msh');
[p,e,t,nVnodes,nPnodes,indices] = convertMeshToSecondOrder(p,e,t);

% Define water dynamic viscosity [Pa*s]
mu = 0.001;

% Generate random permeability field of sand
K = (1e-9)*rand(nVnodes,1);
figure(1)
displaySolution2D(p,t,K,'Porous media permeability');

% Assemble global problem matrix accounting for Darcy law
[M, F] = assembleDiffusionMatrix2D(p,t,K);

% Set high water pressure head at inlet (rho*g*h) with h about 5.5 m
[M,F] = imposeScalarBoundaryCondition2D(p,e,M,F,8,'value',1000*9.81*5.5);
% % Set zero pressure at outlet
[M,F] = imposeScalarBoundaryCondition2D(p,e,M,F,9,'value',0);


% Solve for pressure field
pressure = M\F;

% Display pressure field
figure(2)
displaySolution2D(p,t,pressure,'Pressure field');

% Compute horizontal vx and vertical vy velocity from Darcy's law
gradp = solutionGradient2D(p,t,pressure);
vx = -K/mu.*gradp(:,1);
vy = -K/mu.*gradp(:,2);

% Compute velocity magnitude
vmag = sqrt(vx.^2+vy.^2);

% Display velocity field
figure(3)
displaySolution2D(p,t,vmag,'Velocity magnitude [m/s]');