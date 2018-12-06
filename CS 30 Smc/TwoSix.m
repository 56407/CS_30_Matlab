%TwoSix

diary ('Two_six.out')
a =[2 1;-1 4];
b = [-1 3;0 2];
c = [2;1];
d = eye(2);

%part a)
%result = a + b;
result = a + b;

disp(result);


%part b);
result = a * d;
disp(result);

%part c
result = a * d;
disp(result);

%part d
result = a * c;
disp(result);

%part e
%somethings suppose be illegal
result = a .*c;
disp(result);

%part f
result = a \ b;
disp(result);

%part g
result = a.\b;
disp(result);

%part h
result = a .^ b;
disp(result);


