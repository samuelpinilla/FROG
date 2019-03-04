function x = ifftc(x,varargin)
%IFFTC Centered inverse discrete Fourier transform.
%	IFFTC(X) is the centered inverse discrete Fourier transform of X.
%
%	IFFTC(X,N) is the N-point inverse transform.
% 
%	IFFTC(X,[],DIM) or IFFTC(X,N,DIM) is the inverse discrete Fourier
%	transform of X across the dimension DIM.
%
%	See also FFTC, IFFT, FFT, FFT2, IFFT2, FFTN, IFFTN, FFTSHIFT.

%	v1.0, 6/25/01, Michael Butterfield, <gte881s@prism.gatech.edu>
%	v1.1, 6/26/01, Erik Zeek, <zeekec@mad.scientist.com>
%		Updated help.  Changed to VARARGIN implementation.
%
%	$Revision: 1.2 $ $Date: 2008-08-07 12:56:03 $
%
%	$Log: ifftc.m,v $
%	Revision 1.2  2008-08-07 12:56:03  pam
%	modified to work when the number of points is odd or even
%	
%	Revision 1.1  2006-11-11 00:15:30  pablo
%	CVS server re-installation
%	
%	Revision 1.4  2004/04/01 16:35:04  xg
%	Now accommoates dimensional argument.
%	
%	Revision 1.3  2002/04/18 22:05:29  zeekec
%	Reused the input as the output.
%	
%	Revision 1.2  2001/07/10 01:10:00  zeekec
%	Library cleanup.  Added, deleted, and moved files.
%	
%

if (nargin == 3)
    dim = varargin{2};
    m = size(x);m = m(dim);
    
    if mod(m,2)==1;
        x = fftshift(ifft(ifftshift(x, dim), varargin{:}), dim);
    else
        x = ifftshift(ifft(fftshift(x, dim), varargin{:}), dim);
    end
    
else
    m = size(x);m = m(1); 
    
    if mod(m,2)==1;
        x = fftshift(ifft(ifftshift(x), varargin{:}));
    else
        x = ifftshift(ifft(fftshift(x), varargin{:}));
    end
    
end

    