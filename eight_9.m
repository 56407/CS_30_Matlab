%8.9

t = [0:0.01:10]; %time var

v = 10 * exp((-0.2 + 1i*pi)*t);%func

polar(angle(v),abs(v))

%using polar because of 'pi' modifier
%taking the angle of v because we using a polar function
%take the absolute value of v because its a complex function.