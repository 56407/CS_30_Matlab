%8.14


n = 6; %number of points
a = 2;
b = 1;

t = linspace(-pi,pi,n);
p = linspace(-pi/2,pi/2,n);

[t,p] = meshgrid(t,p); %new meshgrid

x = a * cos(p) * cos(t); %parametric equations
y = b * cos(p) * sin(t);
z = b * sin(p);

plot3(x,y,z); %plot in 3d and surf
surf(x,y,z);
axis equal; % have centered to axis