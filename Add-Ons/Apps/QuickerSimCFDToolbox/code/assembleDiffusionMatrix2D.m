function [D, F] = assembleDiffusionMatrix2D(p,t,nu)

% assembleDiffusionMatrix2D - Assemble global diffusion matrix in 2-D.
%
% This QuickerSim CFD Toolbox function assembles global diffusion matrix.
%
% [D, F] = assembleDiffusionMatrix2D(p, t, nu);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% t - array of finite elements (see help of the importMeshGmsh function for
%     details).
% nu - diffusivity value. This parameter can be entered in one of the
%      following ways:
%      1. As a constant scalar value and then is used in all mesh elements.
%      2. As a 1-by-nelements vector where nelements is the total
%         number of elements in the mesh. In this case nu is assumed
%         to be consant across a single element only and enables
%         implementation of a nonlinear diffusion problem.
%      3. As a nnodes-by-1 vector which specifies diffusivity in each node.
% 
% Output arguments:
% D - assembled, global diffusion matrix with no boundary conditions
%     included.
% F - right-hand side vector of proper size. Because this function does not
%     account for any boundary conditions or source terms, this is always 
%     an all zero vector.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLENAVIERSTOKESMATRIX2D, ASSEMBLESCALARCONVECTIONMATRIX2D,
% ASSEMBLESCALARSOURCETERM2D, ASSEMBLESTOKESMATRIX2D, IMPORTMESHGMSH.

nelements = size(t,2);
nnodes = size(p,2);

if(length(nu)==1)
    nu = repmat(nu,nelements,1);
end

if(length(nu)==nnodes)
    nu = sum(nu(t(1:3,:)))/3;
end

% If first order mesh...
if(size(t,1)==4)
    % Quadrature points
    qPoints = [1/3 1/3];
    weights = [0.5];
    
    nqpoints = size(qPoints,1);
    Nksi = zeros(3,nqpoints);
    Neta = zeros(3,nqpoints);
    
    % Compute derivatives of all linear shape functions in all qPoints
    Nksi(1,:) = -1*ones(1,nqpoints);
    Nksi(2,:) = 1*ones(1,nqpoints);
    Nksi(3,:) = 0*ones(1,nqpoints);
    
    Neta(1,:) = -1*ones(1,nqpoints);
    Neta(2,:) = 0*ones(1,nqpoints);
    Neta(3,:) = 1*ones(1,nqpoints);
    
    % Evaluate nu at centroids of triangles
    
    indRange = 1:9;
    KValues = zeros(9*nelements,1);
    kI = zeros(9*nelements,1);
    kJ = zeros(9*nelements,1);
    
    for el = 1:nelements
        % Compute element level gradient matrices
        elementNodes = t(1:3,el);
        %elementNodes = t(1:6,el);
        x = p(:,elementNodes);
        J = [x(1,2)-x(1,1) x(1,3)-x(1,1);
             x(2,2)-x(2,1) x(2,3)-x(2,1)];

        detJ = det(J);
        L = inv(J);

        Nx = Nksi*L(1,1) + Neta*L(2,1);
        Ny = Nksi*L(1,2) + Neta*L(2,2);

        nuk = repmat(nu(el),1,nqpoints);
        wnuk = weights.*nuk;
        wnuk = repmat(wnuk,3,1);
        %K = nu(el)*(Nx*Nx'+Ny*Ny')*detJ;
        K = (wnuk.*Nx*Nx'+wnuk.*Ny*Ny')*detJ;

        KValues(indRange) = reshape(K,9,1);
        kI(indRange) = repmat(elementNodes,3,1);
        kJ(indRange) = reshape(repmat(elementNodes',3,1),9,1);
        indRange = indRange + 9;
    end
else
    % Quadrature points
    qPoints = [1/6 1/6; 2/3 1/6; 1/6 2/3];
    weights = [1/6 1/6 1/6];

    nqpoints = size(qPoints,1);
    Nksi = zeros(6,nqpoints);
    Neta = zeros(6,nqpoints);

    % Compute derivatives of all quadratic shape functions in all qPoints
    wsp = [1 -3 -3 2 4 2;
                0 -1 0 2 0 0;
                0 0 -1 0 0 2;
                0 4 0 -4 -4 0;
                0 0 0 0 4 0;
                0 0 4 0 -4 -4];

    for i = 1:6
        for k = 1:nqpoints
            ksi = qPoints(k,1);
            eta = qPoints(k,2);

            coords = [0 1 0 2*ksi eta 0;
                      0 0 1 0 ksi 2*eta]';

            Nke = wsp(i,:)*coords;

            Nksi(i,k) = Nke(1);
            Neta(i,k) = Nke(2);
        end
    end

    % Evaluate nu at centroids of triangles

    indRange = 1:36;
    KValues = zeros(36*nelements,1);
    kI = zeros(36*nelements,1);
    kJ = zeros(36*nelements,1);

    for el = 1:nelements
        % Compute element level gradient matrices
        elementPressureNodes = t(1:3,el);
        elementNodes = t(1:6,el);
        x = p(:,elementPressureNodes);
        J = [x(1,2)-x(1,1) x(1,3)-x(1,1);
             x(2,2)-x(2,1) x(2,3)-x(2,1)];

        detJ = det(J);
        L = inv(J);

        Nx = Nksi*L(1,1) + Neta*L(2,1);
        Ny = Nksi*L(1,2) + Neta*L(2,2);

        nuk = repmat(nu(el),1,nqpoints);
        wnuk = weights.*nuk;
        wnuk = repmat(wnuk,6,1);
        %K = nu(el)*(Nx*Nx'+Ny*Ny')*detJ;
        K = (wnuk.*Nx*Nx'+wnuk.*Ny*Ny')*detJ;

        KValues(indRange) = reshape(K,36,1);
        kI(indRange) = repmat(elementNodes,6,1);
        kJ(indRange) = reshape(repmat(elementNodes',6,1),36,1);
        indRange = indRange + 36;
    end
end

D = sparse(kI,kJ,KValues,nnodes,nnodes);

F = zeros(nnodes,1);

end