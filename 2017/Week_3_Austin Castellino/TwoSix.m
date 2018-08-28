%2.6E _ Austin Castellino


diary ('Diary_Two_six.out')
a =[2 1;-1 4];
b = [-1 3;0 2];
c = [2;1];
d = eye(2);

%part a)
%result = a + b;
fprintf('\n Part a) \n');
result = a + b;
disp(result);


%part b);
fprintf('\n Part b) \n');
result = a * d;
disp(result);

%part c
fprintf('\n Part c) \n');
result = a * d;
disp(result);

%part d
fprintf('\n Part d) \n');
result = a * c;
disp(result);

%part e
fprintf('\n Part e) \n');
%somethings suppose be illegal
result = a .*c;
disp(result);

%part f
fprintf('\n Part f) \n');
result = a \ b;
disp(result);

%part g
fprintf('\n Part g) \n');
result = a.\b;
disp(result);

%part h
fprintf('\n Part h) \n');
result = a .^ b;
disp(result);

diary off


