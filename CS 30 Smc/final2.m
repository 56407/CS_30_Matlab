% define x and y axes
x = linspace(0:0.1*pi:10*pi);
y = square(x) * sin(x);
% make plot
plot(x, y, '-.r*');

% set x and y titles
xlabel('x')
ylabel('2x-sin(x)');

plot(x,y,'b--*');