function [ FM, fs ] = create_FM( fc, FM_fm, FM_depth, duration, fs, A, carrier_waveform, FM_waveform)
%[ FM, fs ] = CREATE_FM( fc, FM_fm, FM_depth, duration, fs, A, carrier_waveform, FM_waveform )
% fc : carrier frequency (Hz)
% FM_fm : rate of frequency modulation (Hz)
% FM_depth : beta
% duration : time (in s) (default: 1)
% fs : sampling frequency (Hz) (default: 44100)
% A : amplitude of the signal (default = 1)
% FM_waveform : 'sin' (default), 'triangle', 'square'
% carrier_waveform : 'sin' (default), 'triangle', 'square'
%
% The equation s(t)=1*sin[2*pi*fc*t+beta*sin(2*pi*FM_fm*t+3*pi/2)] describes the
% stimulus. Raised-cosine ramps (20 ms) are applied at the onset and
% at the offset. 
% 
% Leo Varnet 2016 - last modified 2018

%defaults
if nargin<=7
    FM_waveform = 'sin';
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
if nargin<2
    error('Not enough input arguments');
end

if strcmp(carrier_waveform, 'sin')
    function_carrier = @sin;
elseif strcmp(carrier_waveform, 'triangle')
    function_carrier = @triangle;
elseif strcmp(carrier_waveform, 'square')
    function_carrier = @square;
end
if strcmp(FM_waveform, 'sin')
    function_FM = @sin;
elseif strcmp(FM_waveform, 'triangle')
    function_FM = @triangle;
elseif strcmp(FM_waveform, 'square')
    function_FM = @square;
end

% generate t
t = 0:(1/fs):duration;

% stim
stim = function_carrier(2*pi*fc*t+(-1)*FM_depth*function_FM(2*pi*FM_fm*t+3*pi/2));

% ramp
ramp_duration = 0.02;%(1/FM_fm)/2;%duration*0.05;%
t_ramp = 0:(1/fs):ramp_duration;
ramp_up = cos((2*pi*t_ramp(end:-1:1))/(2*max(t_ramp)))/2+0.5;
ramp_down = cos((2*pi*t_ramp)/(2*max(t_ramp)))/2+0.5;
weighting = ones(1,length(t));
weighting(1:length(t_ramp)) = ramp_up;
weighting(end-length(t_ramp)+1:end) = ramp_down;

FM = A*stim.*weighting;

end
