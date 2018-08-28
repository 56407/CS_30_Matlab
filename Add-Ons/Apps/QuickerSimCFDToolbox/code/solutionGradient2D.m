function grad = solutionGradient2D(p, t, u)

% solutionGradient2D - Compute gradient of a scalar field.
%
% This QuickerSim CFD Toolbox function evaluates gradient of a scalar field
% at mesh nodes.
%
% grad = solutionGradient2D(p, t, u)
%
% Input arguments:
% p - array of nodal coordinates described in details in help to the
%     importMeshGmsh function.
% t - array of triangular finite elements described exactly in help to the
%     importMeshGmsh function.
% u - a nnodes-element column vector representing any scalar field (nnodes
%     stands for the total number of mesh nodes). Each element of the u
%     vector represents value of the u scalar field at this particular mesh
%     node.
%
% Output arguments:
% grad - a nnodes-by-2 element matrix whose first column contains values of
%     u field differentiated with respect to x variable and second column
%     with respect to y variable.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: COMPUTEFORCE, COMPUTEMOMENT, EXPORTTOGMSH2D, IMPORTMESHGMSH,
%           SHEARRATE, VORTICITY.

nnodes = size(p,2);
grad = zeros(nnodes,2);
nelements = size(t,2);

if(size(t,1)==4) % first order mesh
    s = zeros(nnodes,1);
    qPoints = [0 0; 1 0; 0 1];
    
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
    
    for el = 1:nelements
        elementNodes = t(1:3,el);
        %elementNodes = t(1:6,el);
        x = p(:,elementNodes);
        J = [x(1,2)-x(1,1) x(1,3)-x(1,1);
             x(2,2)-x(2,1) x(2,3)-x(2,1)];

        %detJ = det(J);
        L = inv(J);

        Nx = Nksi*L(1,1) + Neta*L(2,1);
        Ny = Nksi*L(1,2) + Neta*L(2,2);

        u_local = u(elementNodes);
        ux_local = Nx'*u_local;
        uy_local = Ny'*u_local;

        grad(elementNodes,:) = grad(elementNodes,:) + [ux_local uy_local];
        s(elementNodes) = s(elementNodes) + 1;
    end

    grad = grad./(repmat(s,1,2));
else % second order mesh
    s = zeros(nnodes,1);

    % Nodal locations
    qPoints = [0 0; 1 0; 0 1; 0.5 0; 0.5 0.5; 0 0.5];
    %weights = [1/6 1/6 1/6];

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

    for el = 1:nelements
        elementPressureNodes = t(1:3,el);
        elementNodes = t(1:6,el);
        x = p(:,elementPressureNodes);
        J = [x(1,2)-x(1,1) x(1,3)-x(1,1);
             x(2,2)-x(2,1) x(2,3)-x(2,1)];

        %detJ = det(J);
        L = inv(J);

        Nx = Nksi*L(1,1) + Neta*L(2,1);
        Ny = Nksi*L(1,2) + Neta*L(2,2);

        u_local = u(elementNodes);
        ux_local = Nx'*u_local;
        uy_local = Ny'*u_local;

        grad(elementNodes,:) = grad(elementNodes,:) + [ux_local uy_local];
        s(elementNodes) = s(elementNodes) + 1;
    end

    grad = grad./(repmat(s,1,2));
end
end