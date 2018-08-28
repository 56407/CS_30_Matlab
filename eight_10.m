%8.10


t = (0:0.01:10); %time car

v  = 10*exp((-0.2 + 1i*pi)* t);%func

plot3(real(v), imag(v),t);
%plots the real values of v
%imag(v) helps to show the graphs in a 3d manner
%t is our time




