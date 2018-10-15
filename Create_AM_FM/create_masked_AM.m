function [ AM, fs ] = create_masked_AM( fc, AM_fm, AM_depth, noise_band, sigma_ext, duration, fs )
%[ AM, fs ] = CREATE_MASKED_AM( fc, AM_fm, AM_depth, noise_band, sigma_ext, duration, fs )
% fc : carrier frequency (Hz)
% AM_fm : rate of amplitude modulation (Hz)
% AM_depth : m
% noise_band : low and high cutoff frequencies for the masker
% sigma_ext : standard deviation of the masker
% duration : time (in s) (default: 1)
% fs : sampling frequency (Hz) (default: 44100)
%
% The equation
% s(t)=[1+AM_depth*sin(2*pi*AM_fm*t)+sigma_ext*masker]*sin(2*pi*fc*t)
% describes the stimulus. Raised-cosine ramps (1/2 period) are applied at
% the onset and at the offset. The masker is chosen so that no
% overmodulation occurs.
%
% Leo Varnet & Sarah Attia 2017

%defaults
if nargin<=7
    fs = 44100;
end
if nargin<=6
    duration = 1;
end
if nargin<5
    error('Not enough input arguments');
end

% generate t
t = 0:(1/fs):duration;

% stim
target = AM_depth*sin(2*pi*AM_fm*t);
i = 1;
check = 0;
while i<100 & check==0
    noise = randn(size(t));
    [Bfilt,Afilt] = butter(2, noise_band/(fs/2));
    masker = filter(Bfilt, Afilt, noise);
    masker = masker - mean(masker);
    masker = sigma_ext*masker/sqrt(rms(masker));
    check = ~any(target + masker > 1);
    i=i+1;
end
if i==100
    error('Overmodulation. Try different values for AM_depth and sigma_ext');
end
stim = (1+ target + masker).*sin(2*pi*fc*t);

% ramp
ramp_duration = duration*0.05;
t_ramp = 0:(1/fs):ramp_duration;
ramp_up = cos((2*pi*t_ramp(end:-1:1))/(2*max(t_ramp)))/2+0.5;
ramp_down = cos((2*pi*t_ramp)/(2*max(t_ramp)))/2+0.5;
weighting = ones(1,length(t));
weighting(1:length(t_ramp)) = ramp_up;
weighting(end-length(t_ramp)+1:end) = ramp_down;

AM = stim.*weighting;

end

