function [ gammatones, fc, ERB ] = gammatone_filterbank( fmin, fmax, d, ERBb, fs, n )
%[ gammatones, fc ] = GAMMATONES_FILTERBANK( fmin, fmax, d, ERBb, fs, n )
% Generates a series of gammatones, equally spaced on the ERB scale between
% fmin and fmax with a density of d gammatones per ERB. The ERBb is the
% bandwith of the gammatones, fs the sampling frequency, and n the number
% of points. The filters are stored in matrix gammatones and their center
% frequencies in vector fc.
% if ERBb is set to 0, the approximation of Patterson et al. is used
%
% Leo Varnet 2016

fc = ERBlinspace( fmin, fmax, d );
Nfilters = length(fc);
gammatones = zeros(Nfilters, n);
for i_filter = 1:Nfilters
    clear ERBlim flim fb
    if ERBb>0
        ERBlim = [f2ERB(fc(i_filter))-ERBb/2 f2ERB(fc(i_filter))+ERBb/2];
        flim = ERB2f(ERBlim);
        fb=flim(2)-flim(1);
    else
        fb=0;
    end
    [ gammatones(i_filter,:), ERB(i_filter) ] = irgammatone( fc(i_filter), fb, fs, n );
end
end

