function [idxB, status] = CVodeInitB(fctB, lmmB, nlsB, tB0, yB0, optionsB)
%CVodeInitB allocates and initializes backward memory for CVODES.
%
%   Usage:   IDXB = CVodeInitB ( FCTB, LMMB, NLSB, TB0, YB0 [, OPTIONSB] )
%
%   FCTB     is a function defining the adjoint ODE right-hand side.
%            This function must return a vector containing the current 
%            value of the adjoint ODE righ-hand side.
%   LMMB     is the Linear Multistep Method ('Adams' or 'BDF')
%   NLSB     is the type of nonlinear solver used ('Functional' or 'Newton')
%   TB0      is the final value of t.
%   YB0      is the final condition vector yB(tB0).  
%   OPTIONSB is an (optional) set of integration options, created with
%            the CVodeSetOptions function. 
%
%   CVodeInitB returns the index IDXB associated with this backward
%   problem. This index must be passed as an argument to any subsequent
%   functions related to this backward problem.
%
%   See also: CVodeSetOptions, CVodeInit, CVRhsFnB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/12/05 21:58:18 $

mode = 5;

if nargin < 5
  error('Too few input arguments');
end

if nargin < 6
  optionsB = [];
end

[idxB, status] = cvm(mode, fctB, lmmB, nlsB, tB0, yB0, optionsB);
idxB = idxB+1;
