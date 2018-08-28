function [D11, D12, D22] = assembleCrossDiffusionMatrix2D(p,t,nu)

% assembleCrossDiffusionMatrix2D - This is an internal, undocumented function.

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
    KValues11 = zeros(9*nelements,1);
    KValues12 = zeros(9*nelements,1);
    KValues22 = zeros(9*nelements,1);
    
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
        %K = (wnuk.*Nx*Nx'+wnuk.*Ny*Ny')*detJ;
        
        K11 = (wnuk.*Nx*Nx')*detJ;
        K12 = (wnuk.*Nx*Ny')*detJ;
        K22 = (wnuk.*Ny*Ny')*detJ;

        KValues11(indRange) = reshape(K11,9,1);
        KValues12(indRange) = reshape(K12,9,1);
        KValues22(indRange) = reshape(K22,9,1);
        
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
    KValues11 = zeros(36*nelements,1);
    KValues12 = zeros(36*nelements,1);
    KValues22 = zeros(36*nelements,1);
    
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
        %K = (wnuk.*Nx*Nx'+wnuk.*Ny*Ny')*detJ;

        K11 = (wnuk.*Nx*Nx')*detJ;
        K12 = (wnuk.*Nx*Ny')*detJ;
        K22 = (wnuk.*Ny*Ny')*detJ;

        KValues11(indRange) = reshape(K11,36,1);
        KValues12(indRange) = reshape(K12,36,1);
        KValues22(indRange) = reshape(K22,36,1);
        
        %KValues(indRange) = reshape(K,36,1);
        kI(indRange) = repmat(elementNodes,6,1);
        kJ(indRange) = reshape(repmat(elementNodes',6,1),36,1);
        indRange = indRange + 36;
    end
end

D11 = sparse(kI,kJ,KValues11,nnodes,nnodes);
D12 = sparse(kI,kJ,KValues12,nnodes,nnodes);
D22 = sparse(kI,kJ,KValues22,nnodes,nnodes);

end