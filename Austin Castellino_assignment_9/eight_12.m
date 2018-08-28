%8.12

[x,y] = meshgrid(-1:1/40:1,-2*pi:pi/20:2*pi);
%3d equation
%tricky to read in and program


z = exp(x+1i*y); %our complex function.
figure(1);

mesh(x,y,real(z));
%using real z because we dont want to program the complex points
%of z

