function [elId, b1, b2, b3] = pointInTriangle(p, t, Point)

% pointInTriangle - This is an internal, undocumented function.

x = Point(1);
y = Point(2);

nelements = size(t,2);

elId = 0;

for el = 1:nelements
    node = t(1:3,el);
    
    % Policz odleglosc od prostej (1,2)
    tx = p(1,node(2))-p(1,node(1));
    ty = p(2,node(2))-p(2,node(1));
    L = sqrt(tx^2+ty^2);
    tx = tx/L;
    ty = ty/L;
    nx = -ty;
    ny = tx;
    d1 = (x-p(1,node(1)))*nx + (y-p(2,node(1)))*ny;
    dmax = (p(1,node(3))-p(1,node(1)))*nx + (p(2,node(3))-p(2,node(1)))*ny;
    b3 = d1/dmax;
    
    % Policz odleglosc od prostej (2,3)
    tx = p(1,node(3))-p(1,node(2));
    ty = p(2,node(3))-p(2,node(2));
    L = sqrt(tx^2+ty^2);
    tx = tx/L;
    ty = ty/L;
    nx = -ty;
    ny = tx;
    d2 = (x-p(1,node(2)))*nx + (y-p(2,node(2)))*ny;
    dmax = (p(1,node(1))-p(1,node(2)))*nx + (p(2,node(1))-p(2,node(2)))*ny;
    b1 = d2/dmax;
    
    % Policz odleglosc od prostej (3,1)
    tx = p(1,node(1))-p(1,node(3));
    ty = p(2,node(1))-p(2,node(3));
    L = sqrt(tx^2+ty^2);
    tx = tx/L;
    ty = ty/L;
    nx = -ty;
    ny = tx;
    d3 = (x-p(1,node(3)))*nx + (y-p(2,node(3)))*ny;
    dmax = (p(1,node(2))-p(1,node(3)))*nx + (p(2,node(2))-p(2,node(3)))*ny;
    b2 = d3/dmax;
    
    if(d1>=0 && d2>=0 && d3>=0)
        elId = el;
        %disp([d1 d2 d3 b1 b2 b3]);
        break;
    end
end