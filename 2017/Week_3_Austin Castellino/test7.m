
%function plot3 = test7(x,y,z)
%6.6% and %6.8
%for 6.6 the graph looks like it oscilates between -10 and 10;

% [X,Y,Z] = meshgrid(0:0.1:10);

t = meshgrid(0:0.1:10);
z = (10 * exp(-0.2 + (1j * pi) * t));
plot(z,t);
%surf(v,t);
%6.8% Create a polar plot 
%dependent vars = same as 6.6


polarplot(z,t,'*');

%had a problem graphing with complex numbers 'j'
%made the line into '*' so its easier to show
%im not too sure, if this how the answer was suppose to turn out.




