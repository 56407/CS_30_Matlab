% MATLAB program that displays a tabular form wherein one column represents miles/hr and another column ft/s

% Conversion Values:

% 1 Mile = 5280 feet

% 1 Mile/hr = 5280 / 3600 ft / sec
% 1 Mile/hr = 1.46667 ft/sec

% Input matrix created with values ranging from 0 to 65 with 20 equally spaced values
mph = 0:20:65;

% Creating a two column matrix
% First Column: miles per hour Second Column: Feet per second
for i = 1:length(mph)
% Storing miles per hour value
fps(i,1) = mph(i);
  
% Converting miles per hour to feet per second and storing
fps(i,2) = round(1.4667 * mph(i));

end

% Displaying two column matrix
disp(fps);