function [ BB, AA, fc ] = modulation_filterbank( fmin, fmax, fs, varargin ) % removed N
%[ BB, AA, fc ] = MODULATION_FILTERBANK ( fmin, fmax, fs, NChannels, Qfactor )
%generates a bank of 1st order butterworth filters (default Qfactor=1) between fmin and
%fmax, overlapping at their -3dB point
%
% Leo Varnet and Andrew King 2016

if isempty(varargin) || isempty(varargin{1})
    logfc = log10(fmin):log10((sqrt(5)+1)/(sqrt(5)-1)):log10(fmax);
    Qfactor = 1;
elseif length(varargin)==1
    logfc = linspace(log10(fmin), log10(fmax), varargin{1});
    Qfactor = 1;
elseif length(varargin)==2
    logfc = linspace(log10(fmin), log10(fmax), varargin{1});
    Qfactor = varargin{2};
else
    error('too many input arguments');
end

fc = 10.^(logfc);
for ichan = 1:length(fc)
    flim(ichan,:) = [fc(ichan)*((sqrt(5)-(1/Qfactor))/2) fc(ichan)*((sqrt(5)+(1/Qfactor))/2)];
    [BB(ichan,:),AA(ichan,:)] = butter(2,2*[flim(ichan,:)]/fs);%butter(1,2*[flim(ichan,:)]/fs);
end

end