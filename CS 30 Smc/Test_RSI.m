%rsi formula practice



%its a differential over time/
%period = 14days ,xxHours, or xxMins
%rsi = 100 - (100 / (1+ rs))
%its over +-14 days.
%rs = average gain / average loss
%

rs = input('Please input the rs: ');
%gave me 97.5990
%looking at average loss and gain
%i think i need these two frist ^

rsi =  100 - (100 / (1+ rs));

disp(rsi);

%it doesnt work because other variables have more variables
%Average Gain = [(previous Average Gain) x 13 + current Gain] / 14





