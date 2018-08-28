%meow


function output = timing_7(a1)
vals_of_sort = length(a1);

for ii = 1:vals_of_sort - 1
    ipt = ii;
    
    for jj = ii+1:vals_of_sort
        
        if a1(jj) < a1(ipt)
            ipt = jj;
            
        end
    end
    
    if ii ~= ipt
        temper = a1(ii);
        a1(ii) = a1(ipt);
        a1(ipt) = temper;
    end
end
 output = a1;
 disp(output);