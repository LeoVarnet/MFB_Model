function [ AM, fs, true_envelope ] = create_AM_cloud( fc, fm_range, sigma_ext, duration, fs, A, carrier_waveform, AM_waveform)
%[ AM, fs ] = CREATE_AM_CLOUD( fc, fm, phi, depth, duration, fs, A, carrier_waveform, AM_waveform )
% fc : a vector of N_fc carrier frequencies (Hz), or a vector of edges for
% the N_fc frequency bands in the case where carrier_waveform='noise'
% fm_range : a range [fm_min fm_max] of rates of amplitude modulation (Hz)
% sigma_ext : vector of N_fc depth for the modulations in each frequency channel
% duration : time (in s) (default: 1)
% fs : sampling frequency (Hz) (default: 44100)
% A : amplitude of the signal (default = 1)
% carrier_waveform : 'sin' (default), 'triangle', 'square'
% AM_waveform : 'sin' (default), 'triangle', 'square'
%
% The output is a sum of N_fc tones (carrier frequencies fc), modulated by
% N_fc complex envelopes with fluctuations within the range [fm_min
% fm_max]. Raised-cosine ramps (20 ms) are applied at the onset and at the
% offset.
%
% Leo Varnet 2018

%defaults
if nargin<=7
    AM_waveform = 'sin';
end
if nargin<=6
    carrier_waveform = 'sin';
end
if nargin<=5
    A = 1;
end
if nargin<=4
    fs = 44100;
end
if nargin<=3
    duration = 1;
end
if nargin<=2
    error('Not enough input arguments');
end

if strcmp(carrier_waveform, 'sin')
    function_carrier = @sin;
elseif strcmp(carrier_waveform, 'triangle')
    function_carrier = @triangle;
elseif strcmp(carrier_waveform, 'square')
    function_carrier = @square;
end
if strcmp(AM_waveform, 'sin')
    function_AM = @sin;
elseif strcmp(AM_waveform, 'triangle')
    function_AM = @triangle;
elseif strcmp(AM_waveform, 'square')
    function_AM = @square;
end

if strcmp(carrier_waveform, 'noise')
    N_fc = length(fc)-1;
else
    N_fc = length(fc);
end

% generate t
t = 0:(1/fs):duration;

% generate one envelope for each frequency channel
for i_fc = 1:N_fc
    %%% WARNING: BEATING PROBLEM %%%
    check = 0;
    i=0;
    while i<100 & check==0
        noise = randn(size(t));
        [Bfilt, Afilt] = butter(2, fm_range/(fs/2));
        envelope(i_fc,:) = filter(Bfilt, Afilt, noise);
        envelope(i_fc,:) = envelope(i_fc,:) - mean(envelope(i_fc,:));
        envelope(i_fc,:) = sigma_ext(i_fc)*envelope(i_fc,:)/sqrt(rms(envelope(i_fc,:)));
        check = ~any(envelope(i_fc,:) > 1);
        i=i+1;
    end
    if i==100
        error('Overmodulation. Try different values for AM_depth and sigma_ext');
    end
end

% stim
stim = zeros(1,length(t));

if ~strcmp(carrier_waveform, 'noise')
    for i_fc = 1:N_fc
        %%% WARNING: BEATING PROBLEM %%%
        stim = stim + (1 + envelope(i_fc,:)).*function_carrier(2*pi*fc(i_fc)*t+rand*2*pi);
    end
else
    for i_fc = 1:N_fc
        mod(i_fc,:) = (1 + envelope(i_fc,:)).*randn(size(t));
        [Bfilt,Afilt] = butter(2,fc([i_fc, i_fc+1])/(fs/2));
        mod(i_fc,:) = filtfilt(Bfilt,Afilt,mod(i_fc,:));
        true_envelope(i_fc,:) = abs(hilbert(mod(i_fc,:)));
        stim = stim + mod(i_fc,:);
    end
end

% ramp
ramp_duration = 0.02;%(1/AM_fm)/2;%duration*0.05;
t_ramp = 0:(1/fs):ramp_duration;
ramp_up = cos((2*pi*t_ramp(end:-1:1))/(2*max(t_ramp)))/2+0.5;
ramp_down = cos((2*pi*t_ramp)/(2*max(t_ramp)))/2+0.5;
weighting = ones(1,length(t));
weighting(1:length(t_ramp)) = ramp_up;
weighting(end-length(t_ramp)+1:end) = ramp_down;

AM = A*stim.*weighting;

end
