function F = assembleScalarSourceTerm2D(p,t,sourceValues)

% assembleScalarSourceTerm2D - Compute global source term vector.
%
% This QuickerSim CFD Toolbox assembles global right-hand side vector for a
% given mesh and source term values at specified points.
%
% F = assembleScalarSourceTerm2D(p, t, sourceValues)
%
% Input arguments:
% p - array of nodal coordinates defined in details in help to the
%     importMeshGmsh function.
% t - array of finite elements described in details in help to the
%     importMeshGmsh function.
% sourceValues - a single scalar value of the source term if it is consant
%     in the whole domain or 1-by-n vector, where n can denote either
%     number of elements in the mesh (in this case source value is assumed
%     to be constant across each element) or n can stand for number of
%     nodes in the mesh (in this case source term is interpolated using
%     finite element shape function within the area of each element).
%     sourceValues can also be a handle to the function taking two
%     arguments x and y (in this order) and returning a single scalar value
%     for these given coordinates.
%
% Output arguments:
% F - assembled global right-hand side vector for the given source term.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLEDIFFUSIONMATRIX2D, IMPORTMESHGMSH.

nnodes = size(p,2);
F = zeros(nnodes,1);

% If first order mesh
if(size(t,1)==4)
    qPoints = [1/3 1/3];
    weights = [0.5];
    nqpoints = size(qPoints,1);
    
    % Compute values of all shape functions in all qPoints
    ShapeF = zeros(3,nqpoints);
    
    for i = 1:nqpoints
        ksi = qPoints(i,1);
        eta = qPoints(i,2);
        
        ShapeF(1,i) = 1-ksi-eta;
        ShapeF(2,i) = ksi;
        ShapeF(3,i) = eta;
    end
    
    nelements = size(t,2);
    
    if(isa(sourceValues,'function_handle'))
        % Zrob cos z uchwytem
        sourceAtQuad = zeros(nnodes,1);
        for i = 1:nnodes
            sourceAtQuad(i) = sourceValues(p(1,i),p(2,i));
        end
    else
        % Zrodlo musi byc wektorem wierszowym
        if(size(sourceValues,1)>1)
            sourceValues = sourceValues';
        end

        % Jesli zrodlo jest stale, rozpropaguj je wszystkich punktow
        % kwadartury wszystkich elementow
        if(length(sourceValues)==1)
            sourceAtQuad = repmat(sourceValues,nelements,nqpoints);


        % Jesli zrodlo okreslona dla kazdego elementu, rozpropaguj je
        % wszystkich punktow kwadratury odpowiednich elementow
        elseif(length(sourceValues)==size(t,2))
            sourceAtQuad = repmat(sourceValues',1,nqpoints);


        % Jesli zrodlo okreslono dla kazdego wezla, oblicz dyfuzyjnosci w
        % punktach kwadratury
        elseif(length(sourceValues)==size(p,2))
            Ui = sourceValues(t(1:3,:))';
            sourceAtQuad = Ui*ShapeF;
        end
    end
    
    for el = 1:size(t,2)
        elementPressureNodes = t(1:3,el);
        %elementNodes = t(1:6,el);
        x = p(:,elementPressureNodes);
        J = [x(1,2)-x(1,1) x(1,3)-x(1,1);
             x(2,2)-x(2,1) x(2,3)-x(2,1)];

        detJ = det(J);

        sourcek = sourceAtQuad(el,:);
        fw = sourcek.*weights;
        F(elementPressureNodes) = F(elementPressureNodes) + detJ*ShapeF*fw';
    end
else % second order mesh
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
        sourceAtQuad = zeros(nnodes,1);
        for i = 1:nnodes
            sourceAtQuad(i) = sourceValues(p(1,i),p(2,i));
        end
    else
        % Zrodlo musi byc wektorem wierszowym
        if(size(sourceValues,1)>1)
            sourceValues = sourceValues';
        end

        % Jesli zrodlo jest stale, rozpropaguj je wszystkich punktow
        % kwadartury wszystkich elementow
        if(length(sourceValues)==1)
            sourceAtQuad = repmat(sourceValues,nelements,nqpoints);


        % Jesli zrodlo okreslona dla kazdego elementu, rozpropaguj je
        % wszystkich punktow kwadratury odpowiednich elementow
        elseif(length(sourceValues)==size(t,2))
            sourceAtQuad = repmat(sourceValues',1,nqpoints);


        % Jesli zrodlo okreslono dla kazdego wezla, oblicz dyfuzyjnosci w
        % punktach kwadratury
        elseif(length(sourceValues)==size(p,2))
            Ui = sourceValues(t(1:6,:))';
            sourceAtQuad = Ui*ShapeF;
        end
    end


    for el = 1:size(t,2)
        elementPressureNodes = t(1:3,el);
        elementNodes = t(1:6,el);
        x = p(:,elementPressureNodes);
        J = [x(1,2)-x(1,1) x(1,3)-x(1,1);
             x(2,2)-x(2,1) x(2,3)-x(2,1)];

        detJ = det(J);

        sourcek = sourceAtQuad(el,:);
        fw = sourcek.*weights;
        F(elementNodes) = F(elementNodes) + detJ*ShapeF*fw';
    end
end

end