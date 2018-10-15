function [ AM, fs ] = create_masked_AM_cloud( fc, fm_range, depth, duration, fs, A, carrier_waveform, AM_waveform)
%[ AM, fs ] = CREATE_AM_CLOUD( fc, fm, phi, depth, duration, fs, A, carrier_waveform, AM_waveform )
% fc : a vector of N_fc carrier frequencies (Hz)
% fm_range : a range [fm_min fm_max] of rates of amplitude modulation (Hz)
% depth : modulation index (m)
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

N_fc = length(fc);

% generate t
t = 0:(1/fs):duration;

% stim
stim = zeros(1,length(t));

if ~strcmp(carrier_waveform, 'noise')
    for i_fc = 1:N_fc
        %%% TODO %%%
        check = 0;
        while i<100 & check==0
            noise = randn(size(t));
            [Bfilt,Afilt] = butter(2, fm_range/(fs/2));
            envelope = filter(Bfilt, Afilt, noise);
            envelope = envelope - mean(envelope);
            envelope = sigma_ext*envelope/sqrt(rms(envelope));
            check = ~any(envelope > 1);
            i=i+1;
        end
        if i==100
            error('Overmodulation. Try different values for AM_depth and sigma_ext');
        end
        %stim = (1 + envelope).*sin(2*pi*fc*t);
        stim = stim + (1 + envelope).*function_carrier(2*pi*fc(i_fc)*t+rand*2*pi);
    end
else
    for i_fc = 1:N_fc-1
        envelope = ones(1,length(t));
        for i_fm = 1:N_fm
            for i_phi = 1:N_phi
                envelope = envelope + depth(i_fc, i_fm, i_phi)*function_AM(2*pi*fm(i_fm)*t+phi(i_phi));
            end
        end
        mod = (1+envelope).*randn(size(t));
        [Bfilt,Afilt] = butter(2,fc([i_fc, i_fc+1])/(fs/2));
        mod=filtfilt(Bfilt,Afilt,mod);
        stim = stim + mod;
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
