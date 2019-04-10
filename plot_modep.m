function [ output_args ] = plot_modep( fc, fmc, Emod, varargin )
%PLOT_MODEP plot the excitation pattern of the output of the excitation
%model
%   Use as plot_modep(fc, fmc, Emod, time_EP, decision_statistics) with: 
%   - fc a vector of Nchannels frequency values (spectral filterbank)
%   - fmc a vector of Nmodchannels modulation values (modulation
%     filterbank) 
%   - Emod a Nsamples X Nchannels X Nmodchannels matrix
%   - time_EP the time samples to be considered in the calculation of the
%     excitation pattern (default: all)
%   - decision_statistics the function calculating the statistics to be
%   represented (default: @(x) mean(x.^2,3)) 
%
% Leo Varnet 2016

Nsamples = size(Emod,1);
Nchannels = size(Emod,2);
Nmodchannels = size(Emod,3);
if length(fc)~=Nchannels
    error('fc must be the same size as the first dimension of Emod')
end
if length(fmc)~=Nmodchannels
    error('fc must be the same size as the second dimension of Emod')
end

% defaults
time_EP = 1:Nsamples;
decision_statistics_fun = @(x) mean(x.^2,1); 

if length(varargin)>=2 & ~isempty(varargin{2})
    decision_statistics_fun = varargin{2};
end
if length(varargin)>=1 & ~isempty(varargin{1})
    time_EP = varargin{1};
end

excitation_pattern = squeeze(decision_statistics_fun(Emod(time_EP,:,:)));
figure; imagesc(excitation_pattern');
set(gca,'YDir','normal','XTick', 1:5:Nchannels, 'YTick', 1:Nmodchannels, 'XTickLabels', num2str(floor(fc(1:5:Nchannels)')), 'YTickLabels', num2str(floor(fmc(1:Nmodchannels)')))
%set(h,'EdgeColor','none');
%set(gca,'XScale','log'); set(gca,'YScale','log');%view(2);
%xlim([fc(1) fc(end)]); ylim([fmc(1) fmc(end)]);
xlabel('center frequency (Hz)'); ylabel('modulation frequency (Hz)');
title('modulation excitation pattern');

end

