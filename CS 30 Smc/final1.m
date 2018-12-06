x = {10 8 10 16 15 2 6 7 8 10 12 7 8 9 5 4 2};

initial_fee = 9500;


prompt = '\nwhat is the the price of tution? \n';

i_tution= input(prompt);

prompt = ('what year do you expect to start college?\n');

i_start = input(prompt);


for i = 0;length(x)
    f_cost = (i_tution + (i_tution * (x{i + 1} * .01)));
    
    disp(f_cost);
      
end

% i_start = input(prompt);
% current_fee =0
% current_fee = initial_fee
% total_fee = initial_fee
% for i in range(3):
% current_fee = current_fee + current_fee * increase[i] * .01
% total_fee = total_fee + current_fee
% print "total fee =",total_fee
% 
% Output : total fee = 22974.0