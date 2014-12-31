function status = CVodeReInit(t0, y0, options)
%CVodeReInit reinitializes memory for CVODES
%   where a prior call to CVodeInit has been made with the same
%   problem size N. CVodeReInit performs the same input checking
%   and initializations that CVodeInit does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage: CVodeReInit ( T0, Y0 [, OPTIONS ] ) 
%
%   T0       is the initial value of t.
%   Y0       is the initial condition vector y(t0).  
%   OPTIONS  is an (optional) set of integration options, created with
%            the CVodeSetOptions function. 
%
%   See also: CVodeSetOptions, CVodeInit

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.4 $Date: 2007/12/05 21:58:18 $

mode = 11;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

status = cvm(mode, t0, y0, options);
