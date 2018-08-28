%data set from a previous excerise

X1=[1.1 0.0 -2.1 -3.5 6.0;
0.0 -3.0 -5.6 2.8 4.3 ];



fprintf('\n Data Set: \n ');

size1 = size(X1);

fprintf('No of rows %d: No of columns %d:\n', size1(1), size1(2));
x=[50 60 45 101;90 85 75 102;33 69 89 103];
[row, col] =find(x==max(x(:,2)));
fprintf("Student who got highest score in quiz 2 is with id :%d\n",x(row,4));

fprintf('Average score of each student: \n');
avg=mean(x(:,1:3),2);

%Print average of each student
for i=1:length(avg)
fprintf("Student %d",i);
fprintf(" =%d\n",round(avg(i)));
end