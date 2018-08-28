%x(:,1:3) will contain all the rows with columns from 1 to 3
fprintf('Average score of each student: ');
avg = mean(x(:,1:3),2);

%Print average of each student
for i=1:length(avg)
printf("Student %d",i);
printf(" =%d\n",avg(i));
end

B = 0:1/(2*pi):5*pi;
 disp(B);