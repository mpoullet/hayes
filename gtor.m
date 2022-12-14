function r=gtor(gamma,epsilon)
%GTOR	Inverse Levinson-Durbin recursion.
%----
%USAGE: r=gtor(gamma,epsilon)
%
%	Finds the autocorrelation sequence r(k) from the 
%	reflection coefficients gamma.  If the optional input
%	epsilon is omitted,
%		r=gtor(gamma)
%	then the autocorrelation sequence is normalized so 
%	that r(0)=1.
%
%  see also ATOG, ATOR, GTOA, RTOA, RTOG
%
%---------------------------------------------------------------
% copyright 1996, by M.H. Hayes.  For use with the book 
% "Statistical Digital Signal Processing and Modeling"
% (John Wiley & Sons, 1996).
%---------------------------------------------------------------

   p=length(gamma);
   aa=gamma(1);
   r=[1 -gamma(1)];
   for j=2:p;
      aa=[aa;0]+gamma(j)*[conj(flipud(aa));1];
      r=[r -fliplr(r)*aa];
      end;
   if nargin == 2,
   r = r*epsilon/prod(1-abs(gamma).^2);
   end;
      
