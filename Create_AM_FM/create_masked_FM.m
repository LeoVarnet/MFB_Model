function [ FM, fs ] = create_masked_FM( fc, FM_fm, FM_depth, noise_band, sigma_ext, duration, fs )
%[ FM, fs ] = CREATE_MASKED_FM( fc, FM_fm, FM_depth, duration, fs, A, carrier_waveform, FM_waveform )
% fc : carrier frequency (Hz)
% FM_fm : rate of frequency modulation (Hz)
% FM_depth : beta
% noise_band : low and high cutoff frequencies for the masker
% sigma_ext : standard deviation of the masker
% duration : time (in s) (default: 1)
% fs : sampling frequency (Hz) (default: 44100)
%
% The equation
% s(t)=1*sin[2*pi*fc*t+beta*sin(2*pi*FM_fm*t)+sigma_ext*masker] describes
% the stimulus. Raised-cosine ramps (1/2 period) are applied at the onset
% and at the offset. 
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
target = FM_depth*sin(2*pi*FM_fm*t);

noise = randn(size(t));
[Bfilt,Afilt] = butter(2, noise_band/(fs/2));
masker = filter(Bfilt, Afilt, noise);
masker = masker - mean(masker);
masker = sigma_ext*masker/sqrt(rms(masker));

stim = sin(2*pi*fc*t+target+masker);

% ramp
ramp_duration = duration*0.05;
t_ramp = 0:(1/fs):ramp_duration;
ramp_up = cos((2*pi*t_ramp(end:-1:1))/(2*max(t_ramp)))/2+0.5;
ramp_down = cos((2*pi*t_ramp)/(2*max(t_ramp)))/2+0.5;
weighting = ones(1,length(t));
weighting(1:length(t_ramp)) = ramp_up;
weighting(end-length(t_ramp)+1:end) = ramp_down;

FM = stim.*weighting;

end


