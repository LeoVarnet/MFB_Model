function [out, env] = gammatone_wrapper(S, fs, fc)
%[out, env] = GAMMATONE_WRAPPER(S, fs, fc) filters the sound S (sampling
% frequency: fs) with a 1-ERB-wide gammatone centered on frequency fc.  
% out is the output of the filtering process, env is the envelope of the
% output.
%
% Leo Varnet 2018

C = gammatone_filtering( S, fc, fc, 1, 0, fs, fs/10, 'no', 'no');
out = real(C);

[B,A] = butter(1,2*150/fs);
env = filter(B,A,abs(out));

end