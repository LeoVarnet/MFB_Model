function [ C, E, fc, ERB ] = gammatone_filtering( S, fmin, fmax, d, ERBb, fs, n, display, verbose)
% [ C, E, fc, fb ] = GAMMATONE_FILTERING( S, fmin, fmax, d, ERBb, fs, n, display, verbose )
%   Complex responses C and envelopes E of a gammatone filterbank defined with
%   gammatone_filterbank( fmin, fmax, d, ERBb, fs, n ), applied to sound S.
%   If display='yes', a display of the output is provided.
%   If display='yes', print the processing states.
% %%%% WARNING : ERB doesn't work yet !
%
% Leo Varnet 2016

% default arguments

if nargin<=8
    verbose='no';
end
if nargin<=7
    display='no';
end
if nargin<=6
    n=4410;
end
if nargin<=5
    fs=44100;
end
if nargin<=4
    ERBb=0;
end
if nargin<=3
    d=1;
end
if nargin<=2
    error('Not enough input arguments! Correct syntax for this function is: gammatone_filtering( S, fmin, fmax, ...).');
end

if isyes(verbose)
    fprintf('\ngenerating the gammatone filterbank\n');
end
[gammatones, fc, ERB] = gammatone_filterbank( fmin, fmax, d, ERBb, fs, n );
Nfilters = size(gammatones,1);
Nsamples = length(S);
C = zeros(Nsamples,Nfilters);
for i_filter=1:Nfilters
    if isyes(verbose)
        fprintf(['filtering (gammatone ' num2str(i_filter) ' of ' num2str(Nfilters) ')\n']);
    end
    gammatone_response = filter((gammatones(i_filter,:)),1,S);
    gammatone_response=gammatone_response(:);
    % delay correction
    [~,samples_delay] = max(abs(gammatones(i_filter,:)));
    gammatone_response_correct = [gammatone_response; zeros(samples_delay,1)];
    gammatone_response_correct = gammatone_response_correct(end-Nsamples+1 : end);
    C(:,i_filter) = gammatone_response_correct;
end

E = abs(C);

if isyes(display)
    figure;
    for i_filter=1:Nfilters
        subplot(Nfilters,1,i_filter);
        plot((1:Nsamples)/fs,E(i_filter,:)/max(E(i_filter,:)),'b',(1:Nsamples)/fs,real(C(i_filter,:))/max(real(C(i_filter,:))),'r');
        title(['response of the gammatone filter at ' num2str(fc(i_filter)) ' Hz'])
    end
end
end

