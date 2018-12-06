%8.15

n_p = 20; %number of points
a = 1;
b = 0.5;
tta = linspace(-pi,pi,n_p); %interval of pi on the linspace tta
p = linspace(-pi/2,pi/2,n_p);

[tta,p] = meshgrid(tta,p); %meshgrid

x = a * cos(p) * cos(tta);%parametric equations and functions
y = b * cos(p) * sin(tta);
z = b * sin(p);

surf(x,y,z);

axis equal; %align to the correct axis'

hold on; %since we are ploting 2 functions
%this helps us with the compiler

rad = 2;
tta = linspace(-pi,pi,n_p);
p = linspace(-pi/2,pi/2,n_p);

[tta,p] = meshgrid(tta,p);
x = rad * cos(p) * cos(tta);
y = rad * cos(p) * sin(tta);
z = rad * sin(p);

surf(x,imag(y),z);
alpha(0.5);
axis equal;
hold off;

%the scaling and graph seem to be off.
%i tried using plot3 but it didnt really help me.