% [fid,msg] = fopen('matlabtest.txt','rt');
% if fid > 0
%     [l_x, l_y] = size(matlabtest);%xrow or %yrow
%    disp(l_x);
%    disp(l_y);
%     x = textscan(fid,'%s %s %f %f %f %f %f %f %f',-1); % get the data into x
%     for ii = 1:24 % could use length function to get 24 instead of hard coding
%         total_score = 0;
%         %
%         for jj = 3:9 % same here - use length instead of 9
%             total_score = x{jj}(ii) + total_score; % instead of this loop use sum instead
%             %total_score = sum(ii,jj) + total_score;
%         end
%         x{10}(ii,1) = total_score; % set the final score in a column
%     end
%     status = fclose(fid);
%     average_score = mean(x{10}); % calculate the average of total score
%     y = find(x{10} > average_score); % find returns location of all those above average
%     count = length(y); % returns the number of students above average
%     [n,m] = max(x{10}); % better than sorting, use max to get the top score and its location
%     fprintf('The average score of the class is: %6.2f\n',average_score)
%     fprintf('There are %d students with grades above the average.\n\n',count)
%     fprintf('The highest score belongs to:\n')
%     fprintf('%10s %10s %10.2f\n\n',x{1}{m},x{2}{m},x{10}(m))
%     fprintf('Students who received a grade above the average:\n\n')
%     for ii = 1:count
%         fprintf('%10s %10s %10.2f\n',x{1}{y(ii)},x{2}{y(ii)},x{10}(y(ii)))
%     end
% else
%     disp(msg);
% end

a = (1:10);

total_ax = sum(a,1);
[xrow, ycol] = size(a);
total = cumsum(a);
x{7} (a,1) = total;
%disp(total);
y = max(a);
disp(y);
average_score = mean(x{7});
y = find(x{7} > average_score);
fprintf('The scores above the average: %d\n', y);
y = find(x{7} < average_score);
fprintf('The scores below the average: %d\n', y);
mean(y) = find;

%  for ii = 1:length(a); % could use length function to get 24 instead of hard coding
%         total_score = 0;
%         %
%         for jj = 3:9 % same here - use length instead of 9
%             total_score = x{jj}(ii) + total_score; % instead of this loop use sum instead
%             %total_score = sum(ii,jj) + total_score;
%         end
%         x{7}(ii,1) = total_score; % set the final score in a column
%  end