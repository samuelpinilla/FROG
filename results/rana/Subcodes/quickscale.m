function X = quickscale(Y)
%QUICKSCALE Scales Y
%	QUICKSCALE(Y) Scales the maximum of the absolute value of Y
%	to 1.  Works for both vectors and matricies.
%	
%	X = Y / max(abs(Y(:)));

%	$Revision: 1.1 $ $Date: 2006-11-11 00:15:35 $
%
%	v1.0, 6/27/01, Erik Zeek, <zeekec@mad.scientist.com>
%
%	$Log: quickscale.m,v $
%	Revision 1.1  2006-11-11 00:15:35  pablo
%	CVS server re-installation
%	
%	Revision 1.3  2001/09/21 20:19:39  zeekec
%	Fixed error
%	
%	Revision 1.2  2001/07/10 01:10:00  zeekec
%	Library cleanup.  Added, deleted, and moved files.
%	
%

X = Y / max(abs(Y(:)));