function CG = assembleScalarConvectionMatrix2D(p, t, nu, u, v, supgflag)

% assembleScalarConvectionMatrix2D - Assemble global matrix of a scalar
%                                    convection problem in 2-D.
%
% This QuickerSim CFD Toolbox function assembles global matrix of a scalar
% convection problem. The matrix can be assembled both according to the
% original weak form formulation and with SUPG stabilization.
%
% C = assembleScalarConvectionMatrix2D(p, t, nu, u, v, supgflag)
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% t - array of finite elements (see help of the importMeshGmsh function for
%     details).
% nu - diffusivity of the fluid. For details of specification
%      methods see help for the assembleStokesMatrix2D function. This
%      parameter is needed only in case of stabilization with the SUPG
%      technique. If you assemble nonstabilized form, nu will have no
%      influence on the result - you can use arbitrary value, e.g. 1.
% u - a nVnodes-by-1 vector specifying x-velocity value for each node in
%     the mesh. nVnodes stands for number of velocity nodes in the mesh,
%     i.e. basically all mesh nodes. This vector is used as a convection
%     veloicity for the transported scalar.
% v - a nVnodes-by-1 vector of y-velocity values for each node in the mesh.
%     This vector is used as a convection veloicity for the transported
%     scalar.
% supgflag - a string with one of the following values:
%       'nosupg' - for assembly of convection matrix with no SUPG
%                  stabilization,
%       'supgDoublyAsymptotic' - for assembly of the convection matrix
%                  according to the doubly asymptotic formula for SUPG
%                  stabilization parameter.
%
% Output arguments:
% C - assembled global convection matrix with no boundary conditions
%     imposed.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLEDIFFUSIONMATRIX2D, IMPORTMESHGMSH.

% If first order mesh
if(size(t,1)==4)
    % Quadrature points
    qPoints = [1/3 1/3; 3/5 1/5; 1/5 3/5; 1/5 1/5];
    weights = [-9/32 25/96 25/96 25/96];
    
    nqpoints = size(qPoints,1);
    Nksi = zeros(3,nqpoints);
    Neta = zeros(3,nqpoints);
    ShapeF = zeros(3,nqpoints);
    
    % Compute derivatives of all linear shape functions in all qPoints
    Nksi(1,:) = -1*ones(1,nqpoints);
    Nksi(2,:) = 1*ones(1,nqpoints);
    Nksi(3,:) = 0*ones(1,nqpoints);
    
    Neta(1,:) = -1*ones(1,nqpoints);
    Neta(2,:) = 0*ones(1,nqpoints);
    Neta(3,:) = 1*ones(1,nqpoints);
    
    for k = 1:nqpoints
        ksi = qPoints(k,1);
        eta = qPoints(k,2);

        ShapeF(1,k) = 1-ksi-eta;
        ShapeF(2,k) = ksi;
        ShapeF(3,k) = eta;
    end
    
    nelements = size(t,2);
    nVnodes = max(max(t(1:3,:)));
    indRange = 1:9;
    cValues = zeros(9*nelements,1);
    cI = zeros(9*nelements,1);
    cJ = zeros(9*nelements,1);

    if(size(nu,1)==1 && size(nu,2)==1)
        nu = repmat(nu,nelements,1);
    elseif(size(nu,1)==nVnodes || size(nu,2) == nVnodes)
        % Policz wartosci w srodkach komorek
        nu = sum(nu(t(1:3,:)))/3;
    end
    % Else nu must be equal to the number of elements and you don't need to do
    % anything

    for el = 1:nelements
        % Compute element level gradient matrices
        elementPressureNodes = t(1:3,el);
        elementNodes = t(1:3,el);
        x = p(:,elementPressureNodes);
        J = [x(1,2)-x(1,1) x(1,3)-x(1,1);
             x(2,2)-x(2,1) x(2,3)-x(2,1)];

        detJ = det(J);
        L = inv(J);

        Nx = Nksi*L(1,1) + Neta*L(2,1);
        Ny = Nksi*L(1,2) + Neta*L(2,2);

        U = u(elementNodes);
        V = v(elementNodes);

        uk = ShapeF'*U;
        vk = ShapeF'*V;

        h = Nx.*repmat(uk',3,1) + Ny.*repmat(vk',3,1);
        agradNi = h;
        h = h.*repmat(weights,3,1)*detJ;

        C = ShapeF*h';

        % SUPG stabilisation
        % C = ShapeF*h' + alpha*agradNi*h';
        if(strcmp(supgflag,'supgDoublyAsymptotic'))
            he = sqrt(2*detJ/pi);
            ksi = min(1,he*1/(6*nu(el)));
            alpha = he/(2*1)*ksi;
            C = C + alpha*agradNi*h';
        end

        cValues(indRange) = reshape(C,9,1);
        cI(indRange) = repmat(elementNodes,3,1);
        cJ(indRange) = reshape(repmat(elementNodes',3,1),9,1);
        indRange = indRange + 9;
    end
else % second order mesh
    % Quadrature points
    qPoints = [1/3 1/3; 3/5 1/5; 1/5 3/5; 1/5 1/5];
    weights = [-9/32 25/96 25/96 25/96];

    % if(~strcmp(supgflag,'nosupg'))
    %     qPoints = [0.1012 0.1012; 0.7974 0.1012; 0.1012 0.7974; 0.4701 0.0597; 0.4701 0.4701; 0.0597 0.4701; 0.3333 0.3333];
    %     weights = [0.1259391805448 0.1259391805448 0.1259391805448 0.1323941527885 0.1323941527885 0.1323941527885 0.225]/2;
    %     disp('Przypisalem');
    % end

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

    % Compute values of all quadratic shape functions in all qPoints
    ShapeF = zeros(6,nqpoints);
    for i = 1:6
        for k = 1:nqpoints
            ksi = qPoints(k,1);
            eta = qPoints(k,2);

            coords = [1 ksi eta ksi^2 ksi*eta eta^2]';

            ShapeF(i,k) = wsp(i,:)*coords;
        end
    end

    nelements = size(t,2); 
    nVnodes = max(max(t(1:6,:)));
    indRange = 1:36;
    cValues = zeros(36*nelements,1);
    cI = zeros(36*nelements,1);
    cJ = zeros(36*nelements,1);

    if(size(nu,1)==1 && size(nu,2)==1)
        nu = repmat(nu,nelements,1);
    elseif(size(nu,1)==nVnodes || size(nu,2) == nVnodes)
        % Policz wartosci w srodkach komorek
        nu = sum(nu(t(1:3,:)))/3;
    end
    % Else nu must be equal to the number of elements and you don't need to do
    % anything

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

        U = u(elementNodes);
        V = v(elementNodes);

        uk = ShapeF'*U;
        vk = ShapeF'*V;

        h = Nx.*repmat(uk',6,1) + Ny.*repmat(vk',6,1);
        agradNi = h;
        h = h.*repmat(weights,6,1)*detJ;

        C = ShapeF*h';

        % SUPG stabilisation
        % C = ShapeF*h' + alpha*agradNi*h';
        if(strcmp(supgflag,'supgDoublyAsymptotic'))
            he = sqrt(2*detJ/pi);
            ksi = min(1,he*1/(6*nu(el)));
            alpha = he/(2*1)*ksi;
            C = C + alpha*agradNi*h';
        end

        cValues(indRange) = reshape(C,36,1);
        cI(indRange) = repmat(elementNodes,6,1);
        cJ(indRange) = reshape(repmat(elementNodes',6,1),36,1);
        indRange = indRange + 36;
    end
end

nnodes = size(p,2);
CG = sparse(cI,cJ,cValues,nnodes,nnodes);

end