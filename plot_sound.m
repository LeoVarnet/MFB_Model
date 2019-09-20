function [pdB, f_spec, t_spec] = plot_sound( S, fs, varargin )
%PLOT_SOUND ( S, fs, Freq, Overlap, Nwindow, NFFT, Time, dBlim )
% wavform and spectrographic display of sound S. Spectrogram has parameters
% Freq (frequency range, default [0, 8000]), Overlap (default 0.5), Nwindow
% (default 512), NFFT (default 1024). Time can be used optionnaly to limit
% the sound to an interval. 
% Optionally, you can get the spectrogram as an output:
% [pdB, f_spec, t_spec] = plot_sound( S, fs, varargin )
%
% Leo Varnet - 2018 (last update 2019)

% defaults
dBlim = []; 
Time = []; 
Freq = [0 8000]; 
Overlap = 50/100;
Nwindow = 512;
NFFT = 1024;

if length(varargin)>=6 & ~isempty(varargin{6})
    dBlim = varargin{6};
end
if length(varargin)>=5 & ~isempty(varargin{5})
    Time = varargin{5};
end
if length(varargin)>=4 & ~isempty(varargin{4})
    NFFT = varargin{4};
end
if length(varargin)>=3 & ~isempty(varargin{3})
    Nwindow = varargin{3};
end
if length(varargin)>=2 & ~isempty(varargin{2})
    Overlap = varargin{2};
end
if length(varargin)>=1 & ~isempty(varargin{1})
    Freq = varargin{1};
end

if ~isempty(Time)
    if Time(1)==0
        n1_ech=1;
    else
        n1_ech=Time(1)*fs;
    end
    if Time(2)*fs<length(S)
        n2_ech=Time(2)*fs;
    else
        n2_ech=length(S);
    end
else 
    n1_ech = 1;
    n2_ech = length(S);
end

S=S(n1_ech:n2_ech);
t_S = (1:length(S))/fs;

figure
ht = subplot('Position',[0.1, 0.75, 0.6, 0.2]);
plot(ht, t_S,S);
xlim(ht, [t_S(1) t_S(end)])
ylabel(ht,'Amplitude')
set(ht, 'XTickLabels',{})

load('AuditionColorbar.mat');

[~,f_spec,t_spec,p] = spectrogram(S,Nwindow,Nwindow*Overlap,NFFT,fs);
pdB=10*log10(abs(p));
FreqIndex=find(f_spec>=Freq(1) & f_spec<=Freq(2));

[perio, f_perio] = periodogram(S, 1:length(S), f_spec, fs);
hf = subplot('Position', [0.75, 0.1, 0.2, 0.6]);
perio_dB = 10*log10((perio));%abs
plot(hf, f_perio, perio_dB);
xlim(hf, Freq)
if ~isempty(dBlim)
    ylim(dBlim)
end
ylabel(hf, 'amplitude (dB)');
set(hf, 'XTickLabels',{})
camroll(hf, 90)

hs = subplot('Position',[0.1,0.1, 0.6, 0.6]);
surf(t_spec,f_spec(FreqIndex),pdB(FreqIndex,:),'EdgeColor','none');
surf(hs, t_spec,f_spec(FreqIndex),pdB(FreqIndex,:),'EdgeColor','none');
h=pcolor(t_spec,f_spec(FreqIndex),pdB(FreqIndex,:));

if ~isempty(dBlim)
    caxis(dBlim)
end
set(h,'EdgeColor','none');
axis xy; axis tight; view(0,90);
xlim(hs, [t_S(1) t_S(end)])%xlim(hs,[t_spec(1) t_spec(end)])
ylim(hs, Freq)
xlabel('time (s)'); ylabel('frequency (Hz)');
colormap(Audition);

set(ht, 'Position', [0.1, 0.75, 0.6, 0.2]);
set(hf, 'Position', [0.75, 0.1, 0.2, 0.6]);
set(hs, 'Position', [0.1, 0.1, 0.6, 0.6]);

end

