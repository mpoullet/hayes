function Px = per_smooth(x,win,M,n1,n2)
%SPER                                              [mhh3 5/94]
%       Px = per_smooth(x,n1,n2) computes the smoothed periodogram 
%       of x(n) beginning with x(n1) and ending with x(n2).
%       If n1 and n2 are not specified the periodogram of the 
%	entire sequence will be computed.
%	The window type is specified by win 
%       	1 = rectangular
%		2 = Hamming
%		3 = Hanning
%		4 = Bartlett
%		5 = Blackman
%	and the length in M.
x   = x(:);
if nargin == 3
    n1 = 1;  n2 = length(x);  end;
R  = covar(x(n1:n2),M+1);
r  = [conj(fliplr(R(1,2:M+1))),R(1,1),R(1,2:M+1)];
r'
M  = 2*(M+1)-1;
w  = ones(M,1);
if (win == 2) w = hamming(M);
   elseif (win == 3) w = hanning(M);
   elseif (win == 4) w = bartlett(M);
   elseif (win == 5) w = blackman(M); 
   end;
w
r = r'.*w;
%N  = max(1024,n2-n1+1);
Px = abs(fft(r,1024));
Px(1)=Px(2);
end;