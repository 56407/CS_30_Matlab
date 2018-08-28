%7.3 

%doing hyperbolic values 
%sinh1 cosh1 tanh1

function value = seven_3(str,x)
if nargin < 2
    error('arguements must be 2')
end

if str == 'sinh'
    value = sinh1(x);
end

if str == 'cosh'
    value = cosh1(x);
end

if str == 'tanh'
    value = tanh1(x);
    disp(value);
end

%tanh func
functoin value1 = sinh1(x1)
value1 = (exp(x1) - exp(-x1)) / 2;
end
%sinh func
function value2 = cosh1(x2)
value2 = (exp(x2) + exp(-x2))/2;
end
%cosh
function value3 = tanh1(x3)
value3 = (exp(x3) - exp(-x3)/exp(x3)+exp(-x3));
end








