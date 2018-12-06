% chapter one 1.4E
u = 1;
v = 3;
res = 4*u/(3*v)

% no semicolon on line 4 | just hit enter.
% with semicolon we dont get output
% = 0.4444

%b

res = 2*(v)^-2 / (u + v)^2

% its a little bit tricky with the exponents
% = 1.0385

% c

res = v^3 / ( v^3 - u^3)

% awesome | 1.0385

% D

res =  (4/3)*pi*v^2

%used same variables, u & v for alls
% 37.6991

%1.7E

pwd

%'C:\Users\Austin\Documents\MATLAB\2017'

%1.8E

%mkdir('C:\Users\Austin\Documents\Matlab\')
%File test2 
%Add Folder

%1.9E


('Users\Austin\Documents\MATLAB\2017')


t=-2*pi:pi/10:2*pi;
%calculate | sin(t)|
x=abs(sin(t));
%plot result 
plot (t,x)

diary ChapterOneExcersises.out
%saves output for easy access

diary off
%turn off diary.





