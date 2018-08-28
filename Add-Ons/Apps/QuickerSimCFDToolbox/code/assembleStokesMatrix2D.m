function [SM, F, p, e, t] = assembleStokesMatrix2D(p, e, t, nu, varargin)

% assembleStokesMatrix2D - Assemble global matrix of a stationary Stokes
% problem
%
% This QuickerSim CFD Toolbox function takes the mesh in PET format (see
% help for importMeshGmsh function for detailed description of the format) 
% and kinematic viscocity values and assembles global matrix of the
% corresponding stationary Stokes problem. Any boundary conditions are not
% imposed at this stage of assembly process. They need to be introduced
% later to the matrix with the imposeCfdBoundaryCondition2D function.
%
% [SM, F, p, e, t] = assembleStokesMatrix2D(p, e, t, nu);
% [SM, F, p, e, t] = assembleStokesMatrix2D(p, e, t, nu, sourceTerm);
% 
% Input arguments:
% p, e, t - arrays which store computational mesh (all details of  the PET
%     mesh format described in help for importMeshGmsh function).
% nu - kinematic viscosity of the fluid. This parameter can be entered in
%      one of the following ways:
%      1. As a scalar value and then is used in all mesh elements (this
%         corresponds to the case of laminar flow of a fluid with constant
%         viscosity).
%      2. As a 1-by-nelements vector where nelements is the total
%         number of elements in the mesh. In this case viscosity is assumed
%         to be consant on the area of a single element and enables
%         for example implementation of a nonlinear fluid material model.
% sourceTerm - vector source term in the momentum equation (body force
%      acting on the fluid) which can be passed as a constant vector in the
%      whole domain, constant vector across each element, vector determined
%      in mesh nodes and interpolated automatically on the area of each
%      element with second order shape functions or a function handle. For
%      details of determination consult help of the
%      assembleVectorSourceTerm2D function.
%
% Output argument:
% SM - global Stokes matrix with the following block structure:
%               [ A  B']
%               [ B  0 ]
%      where A has the structure [ D  0 ]
%                                [ 0  D ]
%      with D being diffusion matrix of a single velocity component, B a
%      divergence matrix and transpose of B expresses pressure gradient
%      terms in the momentum equations.
% F -  right hand side vector with 2*nVnodes+nPnodes elements, where
%      nVnodes denotes number of velocity nodes and nPnodes number of
%      pressure nodes. For detailed discussion of velocity and pressure
%      nodes refer to help for the convertToSecondOrderMesh function.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLEVECTORSOURCETERM2D, CONVERTMESHTOSECONDORDER,
%           IMPORTMESHGMSH, IMPOSECFDBOUNDARYCONDITION2D.

% Future description:
%      2. As a vector with length equal to the number of velocity nodes in
%         the mesh. In this case quadratic shape functions are used to
%         interpolate viscosity field inside each of the elements and the
%         diffusion matrix is computed accordingly.


indices = generateIndices2D(p, t);

% Assemble the diffusion matrix
%a = 0;
%f = 0;
%[DM, F] = assempde(model, nu, a, f);
[DM, F] = assembleDiffusionMatrix2D(p,t,nu);

if(max(nu)~=min(nu))
    [D11, D12, D22] = assembleCrossDiffusionMatrix2D(p,t,nu);
end

% Assemble pressure and divergence matrix
[Px, Py] = computeGradientMatrix2D(p, t);

nzmaxDM = nnz(DM);
nzmaxPx = nnz(Px);
nzmaxPy = nnz(Py);
%nzmaxPz = nnz(Pz);

nztotal = 2*nzmaxDM+2*(nzmaxPx+nzmaxPy);

pNodes = t(1:3,:);
vNodesNumber = size(p,2);
pNodesNumber = max(max(pNodes));

if(isempty(varargin))
    F = zeros(2*vNodesNumber+pNodesNumber,1);
    %disp('Varargin pusty')
else
    F = [assembleVectorSourceTerm2D(p,t,varargin{1}); zeros(pNodesNumber,1)];
    %disp('Varargin pelny')
end

% Preallocate
SM = spalloc(2*vNodesNumber+pNodesNumber,2*vNodesNumber+pNodesNumber,nztotal);

if(max(nu)~=min(nu))
    SM(indices.indu,indices.indu) = DM + D11;
    SM(indices.indu,indices.indv) = D12;
    SM(indices.indv,indices.indu) = D12;
    SM(indices.indv,indices.indv) = DM + D22;
else
    SM(indices.indu,indices.indu) = DM;
    SM(indices.indv,indices.indv) = DM;
end

SM(indices.indu,indices.indp) = -Px;
SM(indices.indv,indices.indp) = -Py;

SM(indices.indp,indices.indu) = Px';
SM(indices.indp,indices.indv) = Py';

end