%10.3

function out = sum_any(varargin)

%vars used varagin, nvals, inv, unput, ii-- g_sum, t_sum

out = 0;

%valid inputs
% msg = narginchk(1,Inf, nargin);
% 
% error(msg);


nvals  = size(varargin);

inv_input = 0;

for ii = 1:nvals(2)
    if not (isnumeric(varargin{1,ii}))
       inv_input = 1
       break;
    end
end

g_sum = 0;
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
             g_sum = g_sum + t_sum;
     end
         else 
             disp('Error');
end
         out = g_sum;
     end
     