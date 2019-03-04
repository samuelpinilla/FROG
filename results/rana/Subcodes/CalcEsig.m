function [E]=CalcEsig(P,G)
%ESIG Calculates Esig for FROG.
%	ESIG(P,G) calculates Esig(t,tau) for the given probe pulse, P, 
%	and gate pulse, G.

%	$Revision: 1.1 $ $Date: 2006-11-11 00:15:30 $
%
%	$Log: CalcEsig.m,v $
%	Revision 1.1  2006-11-11 00:15:30  pablo
%	CVS server re-installation
%	
%	Revision 1.3  2002/03/06 23:57:03  zeekec
%	Simplified the calculation.
%	
%	Revision 1.2  2001/10/08 20:38:34  zeekec
%	Removed isrow functions and added temporary variables for the toeplitz arguments.
%	
%	Revision 1.1  2001/08/30 21:53:53  zeekec
%	Renamed Esig to CalcEsig
%	
%	Revision 1.2  2001/07/10 01:10:00  zeekec
%	Library cleanup.  Added, deleted, and moved files.
%
%	v1.2 Michael Butterfield <gte881s@prism.gatech.edu.
%		changed size(E) to length(E) on line 11
%	

error(nargchk(2,2,nargin))

n = length(P);

P = P(:);
G = G(:);

A = [G(n / 2 + 1 : +1 : end); zeros(n / 2 - 0, 1)];
B = [G(n / 2 + 1 : -1 :   1); zeros(n / 2 - 1, 1)];


E = (P * ones(1, n)) .* toeplitz(A, B);
