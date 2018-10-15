function [ fvalues ] = ERBlinspace( fmin, fmax, d )
% [ fvalues ] = ERBLINSPACE( fmin, fmax, d )
%   Generates a row vector of frequencies equally spaced on the ERB scale
%   between fmin and fmax. The number of points is determined by density d
%   (number of points per ERB)
%
% Leo Varnet 2016

ERBmin = f2ERB(fmin);
ERBmax = f2ERB(fmax);
step = 1/d;
ERBvalues = ERBmin:step:ERBmax;
fvalues = ERB2f(ERBvalues);

end

