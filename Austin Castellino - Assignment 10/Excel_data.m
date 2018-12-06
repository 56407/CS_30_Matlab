%excel data 2

%a) all data was imported
%b) Calculate the final score using the formula: (max 10 assignments - i.e. drop the lowest) midterms Final exam Discussion. Add a column to show the final score.


set1 = dataset('File' , 'matlabtest.xls', 'Delimiter',' '); %set1 is set to get a very high index; unusable
%disp(size(matlabtest));%showing that the file and functions work with the use of the of the size func
term_avg = mean(double(set1(0:16:17),2)); % a clever way of writing part b)
%gets the avg, converts to double,'dataset'
%adds the term avg and final exam
sum(term_avg,FinalexamDiscussion);

%c. Calculate the final grade based on 90-100 = A, 80- < 90 B, 70 - < 80 C, 60 - < 70 D, and < 60 F. Add a column to show the final grade.

if FinalScore >= 90
y = 'A';
elseif FinalScore >=80
y = 'B';
elseif FinalScore >=70 
y = 'C';
elseif FinalScore >=60 
y = 'D';
else
y = 'F';
end

%d. Using fprintf, show a table of the final score and grade.
fprintf('%Final Scoref\n',FinalScore);
fprintf('%Gradesf\n',Grade);


%e. Sort the table by the final score, and show the result.
final_result = sortrows(set1,{'FinalScore'},{'ascend'});
% %f.Show the min, max and average final score

stats = grpstats(set1, 'FinalExam', {'min', 'mean', 'max'});