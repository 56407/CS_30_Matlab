function [M, F] = imposeCfdBoundaryCondition2D(p, e, t, M, F, faceIds, bcType, varargin)

% imposeCfdBoundaryCondition2D - Impose boundary conditions for a 2-D flow
%                                problem.
%
% This QuickerSim CFD Toolbox function supplements both global matrix and
% right hand side vector with boundary conditions before solving the
% problem.
%
% [M,F] = imposeCfdBoundaryCondition2D(p,e,t,M,F,faceIds,bcType);
% [M,F] = imposeCfdBoundaryCondition2D(p,e,t,M,F,faceIds,bcType,bcValue);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% e - array of edge elements lying on the boundary of the mesh (see help of
%     the importMeshGmsh function for details).
% t - array of finite elements (see help of the importMeshGmsh function for
%     details).
% M - global Stokes or Navier-Stokes matrix assembled with e.g.
%     assembleStokesMatrix2D or assembleNavierStokesMatrix2D function.
% F - global right hand side vector for the Stokes or Navier-Stokes problem
%     assembled with e.g. assembleStokesMatrix2D or
%     assembleNavierStokesMatrix2D function.
% faceIds - a row vector with ids of the faces (edges in 2D) to which the
%     given boundary condition shall be assigned.
% bcType - A string specifying one of the following BC types:
%       'wall' - this type implies a no-slip flow boundary condition.
%       'slipAlongX' - this type of BC alows x-component of velocity to be
%                      freely computed from the underlying Stokes or
%                      Navier-Stokes equation and the y-velocity is
%                      constrained to be zero at the particular boundary.
%                      This boundary condition type corresponds to a
%                      symmetry BC assumed that the given boundary is
%                      alligned with the x-axis direction.
%       'slipAlongY' - this type of BC alows y-component of velocity to be
%                      freely computed from the underlying Stokes or
%                      Navier-Stokes equation and the x-velocity is
%                      constrained to be zero at the particular boundary.
%                      This boundary condition type corresponds to a
%                      symmetry BC assumed that the given boundary is
%                      alligned with the y-axis direction.
%       'inlet' - this BC type allows for free choice of both x and y
%                 components of velocity vector. Thus, it can be used as a
%                 regular inlet boundary condition, but also as a moving
%                 wall boundary condition. The character of the flow
%                 depends only on how the user defines velocity vector at
%                 the boundary. Velocity vector is always specified as the 
%                 last argument to the function(bcValue) and can be
%                 specified in one of the following ways:
%                 1. As a 2-element row vector with constant values used at
%                 the whole boundary (see Example 2 below).
%                 2. As a handle to the function computing velocity vector
%                 for any x and y point coordinates.
%                 3. As an arbitrary vector field specified for all nodes
%                 in the mesh (however, only appriopriate values from
%                 boundary nodes will be used inside the function). This
%                 way is a straigthforward way in case of very complex flow
%                 fields or vector fields coming from other simulations or
%                 experiments.
%                 For help and examples on ways 2 and 3 go to online
%                 documentation and examples, which you can download freely
%                 at www.quickersim.com/cfd-toolbox-for-matlab/index
%        'rotatingWall' - this BC allows specification of a rotating wall.
%                 The data about position of rotation axis and angular
%                 velocity must be specified in a structure (see example 3
%                 below). This structure is passed to the function as
%                 additional argument.
%        'outflow' - this BC allows to specify a nonzero pressure at
%                 outlet - the average outlet pressure value must be
%                 entered with the bcValue argument.
%     If you do not choose any of the above types for a given boundary, it
%     is automatically assigned an outflow boundary condition (open
%     boundary or do-nothing) boundary condition with zero pressure at
%     outflow.
% bcValue - read above in 'inlet' type of BC.
%
% Output arguments:
% M - original matrix of the flow problem supplemented with new boundary
%     conditions.
% F - original global right hand side vector supplemented with new boundary
%     conditions.
%
% Examples:
%   1. Apply 'wall' boundary condition at egdes 1, 2, 4.
%     [M,F] = imposeCfdBoundaryCondition2D(p,e,t,M,F,[1,2,4],'wall');
%   
%   2. Apply 'inlet' with vx = 1.5 and vy = 0 on the whole edge 3.
%     [M,F] = imposeCfdBoundaryCondition2D(p,e,t,M,F,[3],'inlet',[1.5,0]);
%
%   3. Define rotating wall boundary condition with axis placed at x = 3, y
%      y = 0.5 and angular velocity 30 rad/s (positive values in
%      counter-clockwise direction) at wall with id = 5.
%      rot.axis = [3, 0.5];
%      rot.omega = 30;
%     [M,F] = imposeCfdBoundaryConditino2D(p,e,t,M,F,5,'rotatingWall',rot);
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLENAVIERSTOKESMATRIX2D, ASSEMBLESTOKESMATRIX,
%           IMPORTMESHGMSH, IMPORTSCALARBOUNDARYCONDITION2D.

% Create collection of Dirichlet nodes
nnodes = size(p,2);

% Optimise computational efficiency
edgeNodes = [1 2 4; 2 3 5; 3 1 6];

if(strcmp(bcType,'wall'))
    wallNodes = zeros(1,nnodes);
    
    for faceId = faceIds
        for edge = 1:size(e,2)
            if(e(5,edge)==faceId)
                dnodes = e([1 2 8],edge);
                
                wallNodes(dnodes) = 1;
            end
        end
    end
    
    wallNodes = find(wallNodes);
    
    M(wallNodes,:) = 0;
    M(wallNodes,wallNodes) = eye(length(wallNodes));
    F(wallNodes) = 0;

    M(nnodes+wallNodes,:) = 0;
    M(nnodes+wallNodes,nnodes+wallNodes) = eye(length(wallNodes));
    F(nnodes+wallNodes) = 0;
end

if(strcmp(bcType,'rotatingWall'))
    wallNodes = zeros(1,nnodes);
    
    for faceId = faceIds
        for edge = 1:size(e,2)
            if(e(5,edge)==faceId)
                dnodes = e([1 2 8],edge);
                
                wallNodes(dnodes) = 1;
            end
        end
    end
    
    wallNodes = find(wallNodes);
    mrf = varargin{1};
    
    M(wallNodes,:) = 0;
    M(wallNodes,wallNodes) = eye(length(wallNodes));
    % Calculate vx = -omega*ry (ry = y - axis_y)
    F(wallNodes) = -(mrf.omega)*(p(2,wallNodes)-mrf.axis(2));

    M(nnodes+wallNodes,:) = 0;
    M(nnodes+wallNodes,nnodes+wallNodes) = eye(length(wallNodes));
    % Calculate vy = omega*rx (rx = x - axis_x)
    F(nnodes+wallNodes) = (mrf.omega)*(p(1,wallNodes)-mrf.axis(1));
end

if(strcmp(bcType,'slipAlongX'))
    slipNodes = zeros(1,nnodes);
    
    for faceId = faceIds
        for edge = 1:size(e,2)
            if(e(5,edge)==faceId)
                dnodes = e([1 2 8],edge);
                
                slipNodes(dnodes) = 1;
            end
        end
    end
    
    slipNodes = find(slipNodes);

    M(nnodes+slipNodes,:) = 0;
    M(nnodes+slipNodes,nnodes+slipNodes) = eye(length(slipNodes));
    F(nnodes+slipNodes) = 0;
end

if(strcmp(bcType,'slipAlongY'))
    slipNodes = zeros(1,nnodes);
    
    for faceId = faceIds
        for edge = 1:size(e,2)
            if(e(5,edge)==faceId)
                dnodes = e([1 2 8],edge);
                
                slipNodes(dnodes) = 1;
            end
        end
    end
    
    slipNodes = find(slipNodes);

    M(slipNodes,:) = 0;
    M(slipNodes,slipNodes) = eye(length(slipNodes));
    F(slipNodes) = 0;
end

if(strcmp(bcType,'inlet'))
    inletNodes = zeros(1,nnodes);
    
    for faceId = faceIds
        for edge = 1:size(e,2)
            if(e(5,edge)==faceId)
                dnodes = e([1 2 8],edge);
                inletNodes(dnodes) = 1;
            end
        end
    end
    
    inletNodes = find(inletNodes);
    
    velComp = varargin{1};
    
    if(isa(velComp,'function_handle'))
        % Compute from function
        for i = 1:length(inletNodes)
            v = velComp(p(1,inletNodes(i)), p(2,inletNodes(i)));
            M(inletNodes(i),:) = 0;
            M(inletNodes(i),inletNodes(i)) = 1;
            F(inletNodes(i)) = v(1);

            M(nnodes+inletNodes(i),:) = 0;
            M(nnodes+inletNodes(i),nnodes+inletNodes(i)) = 1;
            F(nnodes+inletNodes(i)) = v(2);
        end
    else
        % Apply BC normally (without function)
        if(size(velComp,1)==1)
            velComp = repmat(velComp,nnodes,1);
        end
        
        M(inletNodes,:) = 0;
        M(inletNodes,inletNodes) = eye(length(inletNodes));
        F(inletNodes) = velComp(inletNodes,1);

        M(nnodes+inletNodes,:) = 0;
        M(nnodes+inletNodes,nnodes+inletNodes) = eye(length(inletNodes));
        F(nnodes+inletNodes) = velComp(inletNodes,2);
    end
end

if(strcmp(bcType,'outflow'))
    % Check pressure value
    pvalue = varargin{1};
    
    for faceId = faceIds
        etmp = e(:,e(5,:)==faceId);
        
        for face = 1:size(etmp,2)
            % Policz normalna
            n1 = etmp(1,face);
            n2 = etmp(2,face);
            n3 = etmp(8,face);
            
            nx = p(2,n2)-p(2,n1);
            ny = p(1,n1)-p(1,n2);
            
            % Policz calke P(t)*n_vec*ShapeF i wpisz do wezlow
            F(n1) = F(n1) - pvalue*nx/6;
            F(n1+nnodes) = F(n1+nnodes) - pvalue*ny/6;
            
            F(n2) = F(n2) - pvalue*nx/6;
            F(n2+nnodes) = F(n2+nnodes) - pvalue*ny/6;
            
            F(n3) = F(n3) - pvalue*nx*2/3;
            F(n3+nnodes) = F(n3+nnodes) - pvalue*ny*2/3;
        end
    end
end

end