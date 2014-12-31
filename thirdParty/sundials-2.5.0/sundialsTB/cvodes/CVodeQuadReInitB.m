function status = CVodeQuadReInitB(idxB, yQB0, optionsB)
%CVodeQuadReInitB reinitializes memory for backward quadrature integration.
%
%   Usage: CVodeQuadReInitB ( IDXB, YS0 [, OPTIONS ] ) 
%
%   IDXB     is the index of the backward problem, returned by
%            CVodeInitB.
%   YQB0     is the final conditions vector yQB(tB0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the CVodeSetQuadOptions function. 
%
%   See also: CVodeSetQuadOptions, CVodeReInitB, CVodeQuadInitB
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/12/05 21:58:18 $

mode = 16;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  optionsB = [];
end
  
idxB = idxB-1;
status = cvm(mode, idxB, yQB0, optionsB);
