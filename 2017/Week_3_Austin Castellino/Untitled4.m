[row, col] =find(x==max(x(:,2)));
fprintf("Student with the highest Quiz 2 score:%d\n",x(row,4));



fprintf('Average score for each Student \n');
avg=mean(x(:,1:3),2);

%Print average of each student
for i=1:length(avg)
fprintf("Student %d",i);

%can do without the round function, however the values are easier to read with it
fprintf(" =%d\n",round(avg(i)));
end