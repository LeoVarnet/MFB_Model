function [ BB, AA, fmc, flim ] = modulation_filterbank( varargin )
%[ BB, AA, fmc, flim ] = MODULATION_FILTERBANK ( fmin, fmax, fs, Qfactor )
%or MODULATION_FILTERBANK ( fmin, fmax, fs, Nmod, Qfactor )
%generates a bank of 1st order butterworth filters between fmin and fmax,
%with an octave spacing (default: overlap at -3db)
%can also be used as 
%[ BB, AA, fmc, flim ] = MODULATION_FILTERBANK ( fmc, fs, Qfactor )
%to generate a filterbank with specific center frequencies fmc
%
% Léo Varnet 2020

% set arguments 
if nargin==3
    fmc = varargin{1};
    fs = varargin{2};
    Qfactor = varargin{3};
end
if nargin==4
    fmin = varargin{1};
    fmax = varargin{2};
    fs = varargin{3};
    Qfactor = varargin{4};
    
    logfmc = log(fmin):log((sqrt(4*Qfactor^2+1)+1)/(sqrt(4*Qfactor^2+1)-1)):log(fmax);
    fmc = exp(logfmc);
end
if nargin==5
    fmin = varargin{1};
    fmax = varargin{2};
    fs = varargin{3};
    Nmod = varargin{4};
    Qfactor = varargin{5};
    
    logfmc = linspace(log(fmin), log(fmax), Nmod);%log(fmin):log((sqrt(4*Qfactor^2+1)+1)/(sqrt(4*Qfactor^2+1)-1)):log(fmax);
    fmc = exp(logfmc);
end

for ichan = 1:length(fmc)
    flim(ichan,:) = fmc(ichan)*sqrt(4+1/Qfactor^2)/2 +  [-1 +1]*fmc(ichan)/Qfactor/2; %sqrt((fmc(ichan)/Qfactor)^2+8*fmc(ichan))/2 + [-1 +1]*(fmc(ichan)/(2*Qfactor));% [fmc(ichan)*((sqrt(5)-(1/Qfactor))/2) fc(ichan)*((sqrt(5)+(1/Qfactor))/2)];
    [BB(ichan,:),AA(ichan,:)] = butter(2,2*[flim(ichan,:)]/fs);%butter(1,2*[flim(ichan,:)]/fs);%
end

% if isempty(varargin) || isempty(varargin{1})
%     logfc = log10(fmin):log10((sqrt(5)+1)/(sqrt(5)-1)):log10(fmax);
%     Qfactor = 1;
% elseif length(varargin)==1
%     logfc = linspace(log10(fmin), log10(fmax), varargin{1});
%     Qfactor = 1;
% elseif length(varargin)==2
%     logfc = linspace(log10(fmin), log10(fmax), varargin{1});
%     Qfactor = varargin{2};
% else
%     error('too many input arguments');
% end
% 
% fc = 10.^(logfc);
% for ichan = 1:length(fc)
%     flim(ichan,:) = [fc(ichan)*((sqrt(5)-(1/Qfactor))/2) fc(ichan)*((sqrt(5)+(1/Qfactor))/2)];
%     [BB(ichan,:),AA(ichan,:)] = butter(2,2*[flim(ichan,:)]/fs);%butter(1,2*[flim(ichan,:)]/fs);
% end

end