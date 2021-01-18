function [ MM, fs ] = create_MM( fc, AM_fm, AM_depth, FM_fm, FM_depth, delta_phi, duration, fs, carrier_waveform, AM_waveform, FM_waveform)
%[ MM, fs ] = CREATE_MM( fc, AM_fm, AM_depth, FM_fm, FM_depth, delta_phi, duration, fs, carrier_waveform, AM_waveform, FM_waveform )
% fc : carrier frequency (Hz)
% AM_fm : rate of amplitude modulation (Hz)
% AM_depth : m
% FM_depth : beta
% FM_fm : rate of frequency modulation (Hz)
% delta_phi : relative phase between AM and FM
% duration : time (in s) (default: 1)
% fs : sampling frequency (Hz) (default: 44100)
% carrier_waveform : 'sin' (default), 'triangle', 'square'
% AM_waveform : 'sin' (default), 'triangle', 'square'
% FM_waveform : 'sin' (default), 'triangle', 'square'
% 
% The equation :
% s(t)=[1+AM_depth*sin(2*pi*AM_fm*t)]*sin[2*pi*fc*t+beta*sin(2*pi*FM_fm*t)]
% describes the stimulus. Raised-cosine ramps (1/2 period) are applied at
% the onset and at the offset.
%
% Leo Varnet 2016

%defaults
if nargin<=11
    FM_waveform = 'sin';
end
if nargin<=10
    AM_waveform = 'sin';
end
if nargin<=9
    carrier_waveform = 'sin';
end
if nargin<=8
    fs = 44100;
end
if nargin<=7
    duration = 1;
end
if nargin<6
    error('Not enough input arguments');
end

if strcmp(AM_waveform, 'sin')
    function_AM = @sin;
elseif strcmp(AM_waveform, 'triangle')
    function_AM = @triangle;
elseif strcmp(AM_waveform, 'square')
    function_AM = @square;
end

% generate t
t = 1/fs:(1/fs):duration;

% stim
FM_stim = create_FM( fc, FM_fm, FM_depth, duration, fs, carrier_waveform, FM_waveform );
MM = (1+AM_depth*function_AM(2*pi*AM_fm*t+delta_phi)).*FM_stim;

end
