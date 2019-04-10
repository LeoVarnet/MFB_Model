function [ E ] = envelope_extraction( s, fs, fc )
% ENVELOPE_EXTRACTION (s, fs, fc)
% extracts the envelopes of s by half-wave rectification and low-pass
% filtering at frequency fc, where s is a Nsamples X Nchan matrix of sounds
% at sampling frequency fs.
%
% Leo Varnet 2016

Nchan = size(s,2);
Nsamples = size(s,1);
E = NaN(Nsamples,Nchan);

[B,A] = butter(1,2*fc/fs);

E = filter(B,A,max(s,0));
    
end