% mfb2.m - applies modulation filterbank				
%
% Usage: [inf, out] = mfb2(in,lmf,umf,den,style,fs)
%
% in  =		input column vector
% lmf =		lowest modulation-filter center frequency,
%      		if 0 the output of a 2nd-order 
%      		Butterworth lowpass filter is additionally computed.
% umf =		highest modulation-filter center frequency,
%      		for typical applications choose umf = 1500.
%      		If lmf = umf only the output of a single filter
%      		at lmf is computed.
% den =		density of overlapping modulation filters,
%      		1 = 3-dB overlap, 0.5 = every second filter, 2 = two times as much filters ... .
% style = 	toggles the overall 150 Hz lowpass on (1) / off (2).
% fs  = 		sampling rate in Hz,
%       		should be greater than 8000 Hz to avoid errors due to the bilinear transform.
%
% inf							= center frequencies of the modulation filters.
% [out1,out2, ...,outn] = each column of martrix out contains the output of
%								  a single modulation filter.
%
% copyright (c) 1999 Stephan Ewert and Torsten Dau, Universitaet Oldenburg

function [inf, out] = mfb2(in,lmf,umf,den,style,fs)

if fs < 8000
   warning('sample rate lower than 8000 Hz')
end

Q = 1;
lpcut = 150;
ex=((3+sqrt(5))/2)^(1/den);

if lmf == 0
   startmf = 4/ex^(den-1);	% dens = 1 or 2
   sw = 1;
elseif lmf > 0
   startmf = lmf;
   sw = 2;
end

if lmf == umf
   mf = startmf;
   sw = 0;
   if lmf == 0
      sw =3;
   end
else
   tmp = fix(log(umf/startmf)/log(ex));
   mf = 0:tmp;
   mf = ex.^mf;
   mf=startmf*mf;
end

b2lpcut = startmf/sqrt(ex);

wlp = 2*pi*lpcut/fs;
wb2lp = 2*pi*b2lpcut/fs;

[b1,a1] = folp(wlp);
[b2,a2] = solp(wb2lp,1/sqrt(2));
%[b2,a2] = folp(wb2lp);

switch style								% switch overall lowpass
   case 1
		outtmp = filter(b1,a1,in);
   case 2
      outtmp = in;
end
   
switch sw
   case 0									% only one modulation filter
		w0 = 2*pi*mf/fs;
   	[b3,a3] = sobp(w0,Q);
      out = filter(b3,a3,outtmp);
      inf = mf;
   case 1									% lowpass and modulation filter(s)
      out = zeros(length(in),length(mf)+1);
		out(:,1) = filter(b2,a2,outtmp);
		for i=1:length(mf)
   		w0 = 2*pi*mf(i)/fs;
   		[b3,a3] = sobp(w0,Q);
   		out(:,i+1) = filter(b3,a3,outtmp);
      end
      inf = [0 mf];
   case 2									% only modulation filters
      out=zeros(length(in),length(mf));
		for i=1:length(mf)
   		w0 = 2*pi*mf(i)/fs;
   		[b3,a3] = sobp(w0,Q);
   		out(:,i) = filter(b3,a3,outtmp);
      end
      inf = mf;
   case 3									% only lowpass
      out = filter(b2,a2,outtmp);
      inf = 0;
end

% subfunctions

% second order bandpass filter
function [b,a] = sobp(w0,Q)

W0 = tan(w0/2);
B0 = W0/Q;

b = [B0; 0; -B0];
a = [1 + B0 + W0^2; 2*W0^2 - 2; 1 - B0 + W0^2];

b = b/a(1);
a = a/a(1);

% first order lowpass filter
function [b,a] = folp(w0);

W0 = tan(w0/2);

b = [W0, W0]/(1 + W0);
a = [1,(W0 - 1)/(1 + W0)];

% second order Butterworth lowpass filter
function [b,a] = solp(w0,Q)

W0 = tan(w0/2);

b = [1; 2; 1];
a = [1 + 1/(Q*W0) + 1/W0^2; -2/W0^2 + 2; 1 - 1/(Q*W0) + 1/W0^2];

b = b/a(1);
a = a/a(1);

% end of mfb2.m







