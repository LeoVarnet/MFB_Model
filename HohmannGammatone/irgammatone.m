function [ k4, ERB, beta, lambda, a, norm, fb ] = irgammatone( fc, fb, fs, n )
% [ k4, beta, lambda, a, norm ] = IRGAMMATONE( fc, fb, fs, n )
%   Impulse response (4th order) of the Gammatone filter of center
%   frequency fc, 3-dB bandwidth fb (see Hohmann, 2002). First n samples of
%   the impulse response. beta is the phase of the filter and lambda its
%   magnitude.
%   if fb is set to 0, the approximation of Patterson et al. is used
%
% Leo Varnet 2016

% phase (beta) of the filter
beta = 2*pi*(fc/fs);

% magnitude (lambda) of the filter
if fb>0
    phi = 2*pi*((fb/2)/fs);
    u = -3/4;
    p = (-2+2*cos(phi)*10^(u/10))/(1-10^(u/10));
    lambda = -p/2-sqrt(((p^2)/4)-1);
else
    a_4 = (pi * 20 * 2 ^ -6);
    c_4 = 2*sqrt(2^(1/4)-1);
    l = 24.7;
    q = 9.265;
    ERB = l + fc/q;
    b = ERB/a_4;
    lambda = exp(-2*pi*b/fs);
    fb = (c_4/a_4)*ERB;
end

% filter coefficient (a) an normalization factor (norm)
a=lambda*exp(1i*beta);
norm = 2*((1-abs(a))^4);

n_ech = (1:n);
k4=a.^(n_ech).*(n_ech.^(3)+6*n_ech.^(2)+11*n_ech+6)./6;
k4=k4*norm;

end
