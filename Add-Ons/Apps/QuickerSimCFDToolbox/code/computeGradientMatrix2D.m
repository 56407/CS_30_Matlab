function [Px, Py] = computeGradientMatrix2D(p, t)

% computeGradientMatrix2D - Compute gradient term and divergence matrices.
%
% This QuickerSim CFD Toolbox function assembles gradient term matrices for
% a Stokes or Navier-Stokes problem and the divergence matrices (which are
% simply transposed matrices of the former ones). The computed matrix is
% assembled according to the formula:
% Px = Integral_over_domain(velocity_shape_f_grad_x * pressure_shape_f)
% Py = Integral_over_domain(velocity_shape_f_grad_y * pressure_shape_f)
% The resulting Px or Py matrix for a single triangular finite element with
% linear shape functions for pressure and second order shape functions for
% velocity is of size 6-by-3.
%
% [Px, Py] = computeGradientMatrix2D(p, t);
%
% Input arguments:
% p - array of node coordinates (detailed description in help for
%     importMeshGmsh function).
% t - array of finite elements (detailed description in help for
%     importMeshGmsh function).
%
% Output arguments:
% Px - assembled global matrix resulting from the discretization of the
%      x-component of pressure gradient in Stokes or Navier-Stokes 
%      equation.
% Py - assembled global matrix resulting from the discretization of the
%      y-component of pressure gradient in Stokes or Navier-Stokes 
%      equation.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLESTOKESMATRIX2D, ASSEMBLENAVIERSTOKESMATRIX2D,
% IMPORTMESHGMSH.

Nksi = zeros(6,3);
Neta = zeros(6,3);
Np = zeros(3,3);

% Quadrature points
qPoints = [1/6 1/6; 2/3 1/6; 1/6 2/3];
weights = [1/3 1/3 1/3];

nqpoints = size(qPoints,1);

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

% Compute values of all linear shape functions in all quadrature points
for k = 1:nqpoints
    ksi = qPoints(k,1);
    eta = qPoints(k,2);

%     Np(1,k) = ksi;
%     Np(2,k) = eta;
%     Np(3,k) = 1-ksi-eta;
    
    Np(1,k) = 1-ksi-eta;
    Np(2,k) = ksi;
    Np(3,k) = eta;
end

% The rest of values depends on the chosen element
nelements = size(t,2);
pxValues = zeros(18*nelements,1);
pyValues = zeros(18*nelements,1);
pxI = zeros(18*nelements,1);
pxyJ = zeros(18*nelements,1);
pyI = zeros(18*nelements,1);

indRange = 1:18;

pNodes = t(1:3,:);
vNodesNumber = size(p,2);
pNodesNumber = max(max(pNodes));

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
    
%     Px = zeros(6,3);
%     Py = zeros(6,3);
%     
%     for i = 1:6
%         for j = 1:3
%             for k = 1:3
%                 Px(i,j) = Px(i,j) + Nx(i,k)*Np(j,k)*detJ/6;
%                 Py(i,j) = Py(i,j) + Ny(i,k)*Np(j,k)*detJ/6;
%             end
%         end
%     end
    
    Px = detJ/6*Nx*Np';
    Py = detJ/6*Ny*Np';
    
    % Insert into sparse matrix
%     indu = indices.indu(elementNodes)';
%     indv = indices.indv(elementNodes)';
%     indw = indices.indw(elementNodes)';
%     indp = indices.indp(elementPressureNodes)';
    
    indu = elementNodes';
    indv = elementNodes';
    indp = elementPressureNodes';
    
    pxValues(indRange) = reshape(Px,18,1);
    pxI(indRange) = repmat(indu',3,1);
    pyI(indRange) = repmat(indv',3,1);
    
    pxyJ(indRange) = reshape(repmat(indp,6,1),18,1);
    pyValues(indRange) = reshape(Py,18,1);
    
    indRange = indRange+18;
end

Px = sparse(pxI,pxyJ,pxValues,vNodesNumber,pNodesNumber);
Py = sparse(pyI,pxyJ,pyValues,vNodesNumber,pNodesNumber);


end