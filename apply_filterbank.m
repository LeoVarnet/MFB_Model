function [ out ] = apply_filterbank( BB, AA, in, filtfilt_option)
%[ out ] = APPLY_FILTERBANK (BB, AA, in)
%   applies the bank of filters described in AA and BB (matrices
%   Nchannels X Ncoef) to the signal in, using filtfilt. in can be a
%   matrix of signals (Nsamples X Nsignals).
%   The output out is a matrix Nsamples X Nsignals X Nchannels
%
% Leo Varnet and Andrew King - 2016 (last update 2019)

if nargin < 4
    filtfilt_option = 'no';
end

Nchannels=size(BB,1);
if size(AA,1)~= Nchannels
    error('Number of lines in AA and BB must be equal');
end
if ndims(in)==1
    in=in(:);
end

Nsamples=size(in,1);
Nsignals=size(in,2);

out=zeros(Nsamples,Nsignals,Nchannels);

for ii=1:Nchannels
    if isyes(filtfilt_option)
        res = filtfilt(BB(ii,:),AA(ii,:),in);
    else
        res = filter(BB(ii,:),AA(ii,:),in);
    end
    out(:,:,ii) = res;
end
end