%7.3
%programming hyperbolic functions

%uses sub functions such as
%sinh1 cosh1 tanh1

function value = seven_3(string,x)

if nargin < 2
    error('arguements must be 2')
end

if string == 'sinh'
    value sinh1(x);
    value = cosh1(x);
end

if string == 'tanh'
    value = tanh1(x);
end

function value1 = sinh1(x1)
value1 = (exp(x1) - exp(-x1))/2;

function value2 = cosh1(x2)
value2 = (exp(x2) + exp(-x2))/2;

function value3 = tanh1(x3)
value3 = (exp(x3) - exp(-x3))/(exp(x3) + exp(-x3);




