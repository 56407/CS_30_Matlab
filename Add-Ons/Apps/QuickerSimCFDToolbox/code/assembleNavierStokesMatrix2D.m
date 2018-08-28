function [NS, F] = assembleNavierStokesMatrix2D(p,e,t,nu,u,v,supgflag,varargin)

% assembleNavierStokesMatrix2D - Assemble global matrix of the stationary
%                                Navier Stokes problem.
%
% This QuickerSim CFD Toolbox function assembles complete global matrix of
% the stationary NavierStokesProblem in 2-D for a known approximation of
% the convection velocity components [u,v]
%
% [NS, F] = assembleNavierStokesMatrix2D(p, e, t, nu, u, v, supgflag);
% [NS, F] = 
%  = assembleNavierStokesMatrix2D(p, e, t, nu, u, v, supgflag, sourceTerm);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% e - array of edge elements lying on the boundary of the mesh (see help of
%     the importMeshGmsh function for details).
% t - array of finite elements (see help of the importMeshGmsh function for
%     details).
% nu - kinematic viscosity of the fluid. For details of specification
%      methods see help for the assembleStokesMatrix2D function.
% u - a nVnodes-by-1 vector specifying x-velocity value for each node in
%     the mesh. nVnodes stands for number of velocity nodes in the mesh,
%     i.e. basically all mesh nodes.
% v - a nVnodes-by-1 vector of y-velocity values for each node in the mesh.
% supgflag - a string with one of the following values:
%       'nosupg' - for assembly of convection matrix with no SUPG
%                  stabilization,
%       'supgDoublyAsymptotic' - for assembly of the convection matrix
%                  according to the doubly asymptotic formula for SUPG
%                  stabilization parameter.
% sourceTerm - vector source term in the momentum equation (body force
%      acting on the fluid) which can be passed as a constant vector in the
%      whole domain, constant vector across each element, vector determined
%      in mesh nodes and interpolated automatically on the area of each
%      element with second order shape functions or a function handle. For
%      details of determination consult hepl of the
%      assembleVectorSourceTerm2D function.
%
% Output arguemnts:
% NS - global, assembled Navier Stokes matrix with no boundary conditions
%      imposed. Boundary conditions to that matrix can be applied with the
%      imposeCfdBoundaryCondition2D function.
% F - right-hand side vector of proper size. Because this function does not
%     account for any boundary conditions or source terms, this is always 
%     an all zero vector.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLEDIFFUSIONMATRIX2D, ASSEMBLESCALARCONVECTIONMATRIX2D,
%           ASSEMBLESCALARSOURCETERM2D, ASSEMBLEVECTORSOURCETERM2D,
%           ASSEMBLESTOKESMATRIX2D, IMPOSECFDBOUNDARYCONDITION2D,
%           IMPORTMESHGMSH.

if(isempty(varargin))
    [SM, F] = assembleStokesMatrix2D(p,e,t,nu);
else
    [SM, F] = assembleStokesMatrix2D(p,e,t,nu,varargin{1});
end

CM = assembleScalarConvectionMatrix2D(p,t,nu,u,v,supgflag);
nztotal = 2*nnz(CM);

indices = generateIndices2D(p, t);

C = spalloc(size(SM,1),size(SM,1),nztotal);
C(indices.indu,indices.indu) = CM;
C(indices.indv,indices.indv) = CM;

NS = SM+C;

end