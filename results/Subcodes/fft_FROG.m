function Esigout=fft_FROG(Esig)
% FFT_FROG Inverse fourier transform Esig.
%   FFT_FROG(Esig) takes in Esig(t,T) and returns the 
%   fourier transform, Esig(w,T).
%
%   See also IFFTC, FFT, IFFT, FFT2, IFFT2, FFTSHIFT.

% v1.0, 6/28/01, Ziyang wang, <gt386x@prism.gatech.edu>
%                Zhenting Dai, <gte859q@prism.gatech.edu>
% v2.0, 7/2/01, Michael Butterfield <gte881s@prism.gatech.edu>
%               Erik Zeek <zeekec@mad.scientist.com>
%       Rewrote FT_Esig.  Vectorized and functionalized.
%
%	$Revision: 1.1 $ $Date: 2006-11-11 00:15:30 $
%
%	$Log: fft_FROG.m,v $
%	Revision 1.1  2006-11-11 00:15:30  pablo
%	CVS server re-installation
%	
%	Revision 1.2  2004/04/01 16:31:04  xg
%	*** empty log message ***
%	
%	Revision 1.1  2001/07/10 01:10:00  zeekec
%	Library cleanup.  Added, deleted, and moved files.
%	
%

Esigout = fftc(Esig, [], 1);