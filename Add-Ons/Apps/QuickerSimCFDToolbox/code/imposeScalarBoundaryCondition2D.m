function [M, F] = imposeScalarBoundaryCondition2D(p, e, M, F, edgeIds, bcType, value, varargin)

% imposeScalarBoundaryCondition2D - Impose Dirichlet, Neumann or Robin
%                   boundary conditions for a 2-D scalar diffusion or
%                   convection-diffusion problem.
%
% This QuickerSim CFD Toolbox function supplements both global matrix and
% right hand side vector with boundary conditions before solving the
% problem.
%
% [M,F] = imposeScalarBoundaryCondition2D(p,e,M,F,edgeIds,bcType,bcValue1);
% [M,F] = imposeScalarBoundaryCondition2D(p,e,M,F,edgeIds,bcType,...
%                                                       bcValue1,bcValue2);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% e - array of edge elements lying on the boundary of the mesh (see help of
%     the importMeshGmsh function for details).
% M - global scalar problem matrix assembled with e.g.
%     assembleDiffusionMatrix2D.
% F - global right hand side vector for the scalar problem assembled with
%     e.g. assembleScalarSourceTerm2D function.
% edgeIds - a row vector with ids of the edges to which the present
%     boundary condition shall be assigned.
% bcType - A string specifying one of the following BC types:
%       'value' - for Dirichlet boundary condition.
%       'flux' - for flux (Neumann type boundary condition).
%       'robin' - for Robin (or third type) boundary condition defined by:
%                 k*grad(u)*normal_vector + bcValue2*u = bcValue1
%                 k is diffusivity used in scalar diffusion equation.
%      For help and examples on using this function and nonconstant
%      boundary conditions refer to Tutorial 5 and 6 in the documentation
%      which you can download freely at:
%      www.quickersim.com/cfd-toolbox-for-matlab/index
%     If you do not choose any of the above types for a given boundary, it
%     is automatically assigned a zero Neumann boundary condition.
% bcValue1 - a single scalar value for Dirichlet or Neumann condition and
%     one of the values in case of Robin condition (see above for details).
% bcValue2 - a single scalar value used only in case of Robin condition.
% Both bcValue1 and bcValue2 can be row vectors with number of elements
% equal to the total number of nodes in the mesh. This allows specifying
% non-constant boundary conditions.
%
% Output arguments:
% M - original matrix of the scalar problem supplemented with new boundary
%     conditions.
% F - original global right hand side vector supplemented with new boundary
%     conditions.
%
% Examples:
%   1. Set temperature of 20 degrees at edges 1, 2, 4.
%     [M,F] = imposeScalarBoundaryCondition2D(p,e,M,F,[1,2,4],'value',20);
%   
%   2. Set heat flux of 50 W/m2 at edge 3.
%     [M,F] = imposeScalarBoundaryCondition2D(p,e,M,F,3,'flux',50);
%
%   3. Impose Robin (convection in case of heat transfer) boundary
%      condition with q = 10*(30-T_wall) at edge 5
% [M,F] = imposeScalarBoundaryCondition2D(p,e,M,F,5,'robin',300,10);
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLENAVIERSTOKESMATRIX2D, ASSEMBLESTOKESMATRIX,
%           IMPORTMESHGMSH, IMPOSECFDBOUNDARYCONDITION2D.

nnodes = size(p,2);
if(~isa(value,'function_handle'))
    if(length(value)==1)
        value = repmat(value,nnodes,1);
    end
else
    % Policz wartosci we wszystkich wezlach, gdzie beda potrzebne
end

if(~isempty(varargin))
    if(~isa(varargin{1},'function_handle'))
        alpha = varargin{1};
        if(length(alpha)==1)
            alpha = repmat(alpha,nnodes,1);
        end
    else
        % Policz wartosci we wszystkich potrzebnych wezlach
    end
end

if(size(e,1) == 7) % 1st order mesh
    if(strcmp(bcType,'value'))
        % Dirichlet condition
        dirichletNodes = extractNodeIdsOnEdges(e,edgeIds);
        M(dirichletNodes,:) = 0;
        M(dirichletNodes,dirichletNodes) = eye(length(dirichletNodes));
        F(dirichletNodes) = value(dirichletNodes);
    elseif(strcmp(bcType,'flux'))
        % Neumann condition
        w = [1 1];
        qpoints = [-sqrt(1/3) sqrt(1/3)];
        % Compute shape function values at quadrature points
        Nik = [0.5*(1-qpoints); 0.5*(1+qpoints)];
        
        for edgeId = edgeIds
            etmp = e(:,e(5,:)==edgeId);
            
            for i = 1:size(etmp,2)
                n1 = etmp(1,i);
                n2 = etmp(2,i);
                Lx = p(1,n2)-p(1,n1);
                Ly = p(2,n2)-p(2,n1);
                L = sqrt(Lx^2+Ly^2);
                J = L/2;
                
                Fk = Nik'*value([n1 n2]);
                
                F(n1) = F(n1) + sum(J*w.*Fk'.*Nik(1,:));
                F(n2) = F(n2) + sum(J*w.*Fk'.*Nik(2,:));
            end
        end
    elseif(strcmp(bcType,'robin'))
        % Robin condition
        w = [1 1];
        qpoints = [-sqrt(1/3) sqrt(1/3)];
        % Compute shape function values at quadrature points
        Nik = [0.5*(1-qpoints); 0.5*(1+qpoints)];
        
        for edgeId = edgeIds
            etmp = e(:,e(5,:)==edgeId);
            
            for i = 1:size(etmp,2)
                n1 = etmp(1,i);
                n2 = etmp(2,i);
                Lx = p(1,n2)-p(1,n1);
                Ly = p(2,n2)-p(2,n1);
                L = sqrt(Lx^2+Ly^2);
                J = L/2;
                
                Fk = Nik'*value([n1 n2]);
                alphak = Nik'*alpha([n1 n2]);
                
                F(n1) = F(n1) + sum(J*w.*Fk'.*Nik(1,:));
                F(n2) = F(n2) + sum(J*w.*Fk'.*Nik(2,:));
                
                % Aij = integral( alpha * J * Ni * Nj )
                A11 = sum(J*w.*alphak'.*Nik(1,:).*Nik(1,:));
                A12 = sum(J*w.*alphak'.*Nik(1,:).*Nik(2,:));
                A21 = A12;
                A22 = sum(J*w.*alphak'.*Nik(2,:).*Nik(2,:));
                
                % Add contributions to the matrix
                M(n1,n1) = M(n1,n1) + A11;
                M(n1,n2) = M(n1,n2) + A12;
                M(n2,n1) = M(n2,n1) + A21;
                M(n2,n2) = M(n2,n2) + A22;
            end
        end
    end
else
    % Second order mesh
    if(strcmp(bcType,'value'))
        % Dirichlet condition
        dirichletNodes = extractNodeIdsOnEdges(e,edgeIds);
        M(dirichletNodes,:) = 0;
        M(dirichletNodes,dirichletNodes) = eye(length(dirichletNodes));
        F(dirichletNodes) = value(dirichletNodes);
    elseif(strcmp(bcType,'flux'))
        % Neumann condition
        w = [5/9 8/9 5/9];
        qpoints = [-sqrt(0.6) 0 sqrt(0.6)];
        % Compute shape function values at quadrature points
        Nik = [-0.5*qpoints.*(1-qpoints); 0.5*qpoints.*(1+qpoints); 1-qpoints.^2];
        
        for edgeId = edgeIds
            etmp = e(:,e(5,:)==edgeId);
            
            for i = 1:size(etmp,2)
                n1 = etmp(1,i);
                n2 = etmp(2,i);
                n3 = etmp(8,i);
                Lx = p(1,n2)-p(1,n1);
                Ly = p(2,n2)-p(2,n1);
                L = sqrt(Lx^2+Ly^2);
                J = L/2;
                
                Fk = Nik'*value([n1 n2 n3]);
                
                F(n1) = F(n1) + sum(J*w.*Fk'.*Nik(1,:));
                F(n2) = F(n2) + sum(J*w.*Fk'.*Nik(2,:));
                F(n3) = F(n3) + sum(J*w.*Fk'.*Nik(3,:));
            end
        end
    elseif(strcmp(bcType,'robin'))
        % Robin condition
        w = [5/9 8/9 5/9];
        qpoints = [-sqrt(0.6) 0 sqrt(0.6)];
        % Compute shape function values at quadrature points
        Nik = [-0.5*qpoints.*(1-qpoints); 0.5*qpoints.*(1+qpoints); 1-qpoints.^2];
        
        for edgeId = edgeIds
            etmp = e(:,e(5,:)==edgeId);
            
            for i = 1:size(etmp,2)
                n1 = etmp(1,i);
                n2 = etmp(2,i);
                n3 = etmp(8,i);
                Lx = p(1,n2)-p(1,n1);
                Ly = p(2,n2)-p(2,n1);
                L = sqrt(Lx^2+Ly^2);
                J = L/2;
                
                Fk = Nik'*value([n1 n2 n3]);
                alphak = Nik'*alpha([n1 n2 n3]);
                
                F(n1) = F(n1) + sum(J*w.*Fk'.*Nik(1,:));
                F(n2) = F(n2) + sum(J*w.*Fk'.*Nik(2,:));
                F(n3) = F(n3) + sum(J*w.*Fk'.*Nik(3,:));
                
                % Aij = integral( alpha * J * Ni * Nj )
                A11 = sum(J*w.*alphak'.*Nik(1,:).*Nik(1,:));
                A12 = sum(J*w.*alphak'.*Nik(1,:).*Nik(2,:));
                A13 = sum(J*w.*alphak'.*Nik(1,:).*Nik(3,:));
                A21 = A12;
                A22 = sum(J*w.*alphak'.*Nik(2,:).*Nik(2,:));
                A23 = sum(J*w.*alphak'.*Nik(2,:).*Nik(3,:));
                A31 = A13;
                A32 = A23;
                A33 = sum(J*w.*alphak'.*Nik(3,:).*Nik(3,:));
                
                % Add contributions to the matrix
                M(n1,n1) = M(n1,n1) + A11;
                M(n1,n2) = M(n1,n2) + A12;
                M(n1,n3) = M(n1,n3) + A13;
                M(n2,n1) = M(n2,n1) + A21;
                M(n2,n2) = M(n2,n2) + A22;
                M(n2,n3) = M(n2,n3) + A23;
                M(n3,n1) = M(n3,n1) + A31;
                M(n3,n2) = M(n3,n2) + A32;
                M(n3,n3) = M(n3,n3) + A33;
            end
        end
    end
end

end