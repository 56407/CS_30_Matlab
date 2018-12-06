%10.3

function out = ten_3(varargin)


out = 0; %intialize the function by giving it the value of 0


nvals  = size(varargin); %how many arguements are allowed?

inv_input = 0;

for ii = 1:nvals(2)
    if not (isnumeric(varargin{1,ii})) %checking for valid arguements, 'isnumeric'
       inv_input = 1
       break;
    end
end

sum = 0; %intialize
 %start adding stuff
 
 if inv_input == 0
     for ii = 1:nvals(2) 
         t_sum = 0;
         [m,n] = size(varargin{1,ii});
         if n == 0
             t_sum = sum(varargin{1, ii});
         elseif n >= 1
                 t_sum = sum(sum(varargin{1,ii}));
         end
             sum = sum + t_sum;
     end
         else 
             disp('Error');
end
         out = sum;
end
     
%pretty confusing
%definetly need to take another look at the adding methods
%check how and what 'varargin' does
     