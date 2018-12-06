%10.6 

function out = ten_six(p1,p2)

out = 0;

if isnumeric(p1.x) && isnumeric(p1.y)&& isnumeric(p2.x)&&isnumeric(p2.y)
    out = sqrt((p2.x - p1.x)^2 + (p2.y - p1.y)^2);
end
end


%straight forward
%made a function that calculates the distance
%for the function it is very important to intialize the output
%using the isnumeric modifier to check the arguements of the func
%since we are using p1 and p2 in the code those vars are already occuppied
%and can't be used when writing and using the function

%tricky syntax and interesting that we had to explicitly - 
%make x and y coordinates for both our points, p1 & p2;

