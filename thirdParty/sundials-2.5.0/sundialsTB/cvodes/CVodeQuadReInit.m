function status = CVodeQuadReInit(yQ0, options)
%CVodeQuadReInit reinitializes CVODES's quadrature-related memory
%   assuming it has already been allocated in prior calls to CVodeInit 
%   and CVodeQuadInit.
%
%   Usage: CVodeQuadReInit ( YQ0 [, OPTIONS ] ) 
%
%   YQ0      Initial conditions for quadrature variables yQ(t0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the CVodeSetQuadOptions function. 
%
%   See also: CVodeSetQuadOptions, CVodeQuadInit

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/12/05 21:58:18 $

mode = 12;

if nargin < 1
  error('Too few input arguments');
end

if nargin < 2
  options = [];
end
  
status = cvm(mode, yQ0, options);
