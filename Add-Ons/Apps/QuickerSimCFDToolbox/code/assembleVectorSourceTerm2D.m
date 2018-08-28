function F = assembleVectorSourceTerm2D(p,t,sourceValues)

% assembleVectorSourceTerm2D - Compute global source term vector.
%
% This QuickerSim CFD Toolbox assembles global right-hand side vector for a
% given mesh and source term values at specified points.
%
% F = assembleVectorSourceTerm2D(p, t, sourceValues)
%
% Input arguments:
% p - array of nodal coordinates defined in details in help to the
%     importMeshGmsh function.
% t - array of finite elements described in details in help to the
%     importMeshGmsh function.
% sourceValues - a two-element row vector of the source term components if 
%     it is consant in the whole domain or a 2-by-n vector, where n can
%     denote either number of elements in the mesh (in this case source
%     value is assumed to be constant across each element, the first row
%     determines x-components of the vector source term on each triangle
%     and second row y-components at each of the triangles) or n can stand
%     for the number of nodes in the mesh (in this case source term is
%     interpolated using shape functions within the area of
%     each element and analogically row 1 contains x-components and row 2
%     y-components of the vector source term). Also a handle to the
%     function may be passed as the sourceValues argument. The function
%     should accept two scalars x and y as arguments and return a 2-element
%     row vector as resulting value. For details visit
%     www.quickersim.com/cfd-toolbox-for-matlab/index and download
%     additional free of charge examples and tutorials.
%
% Output arguments:
% F - assembled global right-hand side vector for the given source term.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLEDIFFUSIONMATRIX2D, ASSEMBLESCALARSOURCETERM2D,
%           IMPORTMESHGMSH.

nnodes = size(p,2);
F = zeros(2*nnodes,1);

% Quadrature points
qPoints = [1/6 1/6; 2/3 1/6; 1/6 2/3];
weights = [1/6 1/6 1/6];

nqpoints = size(qPoints,1);

% Compute values of all quadratic shape functions in all qPoints
wsp = [1 -3 -3 2 4 2;
            0 -1 0 2 0 0;
            0 0 -1 0 0 2;
            0 4 0 -4 -4 0;
            0 0 0 0 4 0;
            0 0 4 0 -4 -4];
        
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

if(isa(sourceValues,'function_handle'))
    % Zrob cos z uchwytem
    sourceAtQuadX = zeros(nnodes,1);
    sourceAtQuadY = zeros(nnodes,1);
    for i = 1:nnodes
        v = sourceValues(p(1,i),p(2,i));
        sourceAtQuadX(i) = v(1);
        sourceAtQuadY(i) = v(2);
    end
else
    % Zrodlo musi byc wektorem wierszowym
%     if(size(sourceValues,1)>1)
%         sourceValues = sourceValues';
%     end

    % Jesli zrodlo jest stale, rozpropaguj je wszystkich punktow
    % kwadartury wszystkich elementow
    if(size(sourceValues,1)==1)
        %sourceValues(1)
        sourceAtQuadX = repmat(sourceValues(1),nelements,nqpoints);
        sourceAtQuadY = repmat(sourceValues(2),nelements,nqpoints);


    % Jesli zrodlo okreslona dla kazdego elementu, rozpropaguj je
    % wszystkich punktow kwadratury odpowiednich elementow
    elseif(size(sourceValues,2)==size(t,2))
        sourceAtQuadX = repmat(sourceValues(1,:)',1,nqpoints);
        sourceAtQuadY = repmat(sourceValues(2,:)',1,nqpoints);

    % Jesli zrodlo okreslono dla kazdego wezla, oblicz dyfuzyjnosci w
    % punktach kwadratury
    elseif(size(sourceValues,2)==size(p,2))
        sx = sourceValues(1,:);
        sy = sourceValues(2,:);
        Uix = sx(t(1:6,:))';
        Uiy = sy(t(1:6,:))';
%         disp(size(Uix))
%         disp(size(ShapeF))
        sourceAtQuadX = Uix*ShapeF;
        sourceAtQuadY = Uiy*ShapeF;
    end
end


for el = 1:size(t,2)
    elementPressureNodes = t(1:3,el);
    elementNodes = t(1:6,el);
    x = p(:,elementPressureNodes);
    J = [x(1,2)-x(1,1) x(1,3)-x(1,1);
         x(2,2)-x(2,1) x(2,3)-x(2,1)];

    detJ = det(J);
    
    sourcekx = sourceAtQuadX(el,:);
    fwx = sourcekx.*weights;
    sourceky = sourceAtQuadY(el,:);
    fwy = sourceky.*weights;
    F(elementNodes) = F(elementNodes) + detJ*ShapeF*fwx';
    F(elementNodes+nnodes) = F(elementNodes+nnodes) + detJ*ShapeF*fwy';
end

end