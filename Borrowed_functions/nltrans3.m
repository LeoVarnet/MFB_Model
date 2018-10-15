% nltrans3.m - nonlinear transducer -
%
% Usage: out = nltrans(in,knee,exp)
%
% in = input time domain signal
% knee = compression knee point, e.g., 0.1
% exp = compression in dB/dB, e.g., 0.25
% smooth = ammount of smoothing around knee point (0 = no smoothing)
%
% out = output waveform

% based on my nltrans.m (1999) and nltrans2.m (2007)
% Stephan Ewert, 1999-2016

function out = nltrans3(in,knee,exp,smooth)


if smooth == 0

    idx = find(abs(in) > knee);
    out = in;
    out(idx) = sign(out(idx)).*(abs(out(idx)).^exp / knee^exp * knee);

else
    
    sExp = 1/smooth;
    % root of sum of squares smoothing of the inverted function (then
    % inverted back) 
    % FIXME: hacky. it has division by zero bit Matlab does not complain
    out =sign(in).*1./( (1./abs(in)).^sExp + (1./( abs(in).^exp / knee^exp * knee  )).^sExp ).^(1/sExp);

end

% eof