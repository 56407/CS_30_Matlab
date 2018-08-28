a = (1:10);

total_ax = sum(a,1);
[xrow, ycol] = size(a);
total = cumsum(a);
x{7} (a,1) = total;
%disp(total);
y = max(a);
fprintf('The max of the matrix is: %d \n', y);

y = min(a);

fprintf('The min of the matrix is: %d\n', y);
average_score = mean(x{7});
y = find(x{7} > average_score);
fprintf('The scores above the average: %d\n', y);
y = find(x{7} < average_score);
fprintf('The scores below the average: %d\n', y);
% 
% t_average_score = mean(average_score(a));

t_average_score = mean(a);
disp(t_average_score);