% Week 3 homework
diary ('week_3_data.out')

%Data set from 2.4E
array1=[1.1 0.0 -2.1 -3.5 6.0;
    0.0 -3.0 -5.6 2.8 4.3;
    2.1 0.3 0.1 -0.4 1.3;
    -1.4 5.1 0.0 1.1 -3.0];

fprintf('\n Data: \n');
disp(array1());
% display the contents of the array;
%row one column one
%r1, c1
%disp(array1([1 2],:))

%find the max of the array
%2d array, show the max and the postion it is 
%in the array
fprintf('\n     Show the max values of the array: \n');
max(array1);
[a,b] = max(array1)
fprintf('c shows the position of the max: \n');
% no semi colon, limited scope*
%the array

%min
fprintf('\n     Show the min values of the array: \n');
min(array1);
[c,d] = min(array1)
fprintf('d shows the position of the min: \n');

%average | mean
fprintf('\n     Show the average/mean values of the array: \n');
M_a1 = mean(array1)
disp(M_a1);

