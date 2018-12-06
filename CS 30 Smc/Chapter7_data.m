%chapter 7 data set

  
   a1 = [4.91 8.18;
         3.84 7.49;
         2.41 7.11;
         2.62 6.15;
         3.78 6.62;
         0.52 3.30;
         1.83 2.05;
         2.01 2.83;
         0.28 1.16;
         1.08 0.52;
         0.94 0.21;
         0.59 1.73;
         0.69 3.96;
         3.04 4.26;
         1.01 6.75;
         3.60 6.67;
         4.53 7.70;
         6.13 7.31;
         4.43 9.05;
         6.13 7.31;
         4.43 9.05;
         4.12 10.95];
     
         
   %disp(a1);
   %disp(a1); shows us the columns and rows to 1 - 20
   
temp_a1 = 0;
temp_a2 = 0;
tic
    disp(sort(a1, 1));
    fprintf('sort function: ');
temp_a1 = toc;
disp(temp_a1); % time it took for this sort

tic
fprintf('Max and Min function time: ');
    disp(max(a1));
    disp(min(a1));
    temp_a2 = toc;
    disp(temp_a2); %displays time it took to get the max and min vals of the array
    
    if temp_a1 > temp_a2
        fprintf("\nThe sorting method is faster");
    else 
        fprintf("\nThe min and max functions are faster");
    end
        
        
   %this is as far as I got :(
   %I understand that the tic - toc  func is used to time the 
   %sorted functions.
   %I managed to get the tic toc function to work, However 
   %I'm unsure if I answered the question correctly.
        
        
% misc code\\\\\\\\\\\\\\\\\\\\ignore
%    %plot(a1);
%    function output = output(a1);
%  output = a1;
%  temp_a1 = a1;
%    
%    
%    temp_a1 = initial_data;
%    
%    tic;
%    output = sort(a1);
%    time_passed_a1  = toc;
%    temp_a1 = initial_data;
%    
%    tic;
%    output = sort(temp_a1);
%    time_passed_a2 = toc;
%    
%      disp(num2str(time_passed_a1);
%      disp(num2str(time_passed_a2);

% 
% for l = 1:1:length(a1)
% a1 =   [ 4.91 3.84 2.41 2.62 3.78 0.52 1.83 2.01 0.28 1.08 0.94 0.59...
%        0.69 3.04 1.01 3.60 4.53 6.13 4.43 4.12]; 
% a2 = [8.18 7.49 7.11 6.15 6.62 3.30 2.05 2.83 1.16 0.52 0.21 1.73 3.96...
%        4.26 6.75 6.67 7.70 7.31 9.05 10.95];
%    fprintf(l,(" : "), a1(l), a2(l));
%   
%     %disp(l);
%     
% end
% a1 =   [ 4.91 3.84 2.41 2.62 3.78 0.52 1.83 2.01 0.28 1.08 0.94 0.59...
%        0.69 3.04 1.01 3.60 4.53 6.13 4.43 4.12;
%        
%        8.18 7.49 7.11 6.15 6.62 3.30 2.05 2.83 1.16 0.52 0.21 1.73 3.96...
%        4.26 6.75 6.67 7.70 7.31 9.05 10.95];
   %a1_t = array2table(a1,'x', 'y');
     
   
   
   
 
   
   
   
   
   
   