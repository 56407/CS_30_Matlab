%assignment 11
%From this file, and using the assignment specs from last week, show the first and last name of the student with the highest score.
%Show the list of students (first and last name) and the count of those students who's final score is above the average.
%taking this one a bit slower than the others


%disp(matlabtest); %confirmed that the file is present and working
%upon running, data is formatted and properly organized.

%now lets see if we can sort a column
%disp(double(mean(matlabtest(:,1:2))));
%a_t = table('Last', 'Frist', 'a1', 'a2', 'a3', 'a4' ,'a5' , 'FinalScore');
%disp(a_t);
a = sortrows(matlabtest(:, :)); %the dataset works!!!
a = sortrows(a,'Chamkit','ascend'); %sorting by last name
disp(a);
%a(:,'VarName3') = []; %how to delete a column
%a(:,'VarName2') = [V2]; %how to rename a column.
%or of similar methods
fprintf('\n set of a, sorted by name \n');
disp(a);
fprintf('\n set of a2 \n');


%disp(a); %checking to see that everything is fine

%now working with a2, where we are sorting the scores column
a2 = sortrows(a(1:24, 1:8)); %using (7:9) shows 3 sorted rows
%note that the subscript can be confusing because we are working with a
%table and not a like a double or int, i.e standard array
%review
%also with the subscript, it becomes easier to work with values


fprintf('\n Final Scores \n');
disp(a2); %finalscores % Scores listed



