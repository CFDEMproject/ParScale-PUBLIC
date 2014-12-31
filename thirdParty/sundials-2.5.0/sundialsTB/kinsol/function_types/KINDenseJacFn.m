%KINDenseJacFn - type for user provided dense Jacobian function.
%
%   The function DJACFUN must be defined as 
%        FUNCTION [J, FLAG] = DJACFUN(Y,FY)
%   and must return a matrix J corresponding to the Jacobian of f(y).
%   The input argument FY contains the current value of f(y).
%   If a user data structure DATA was specified in KINInit, then
%   DJACFUN must be defined as
%        FUNCTION [J, FLAG, NEW_DATA] = DJACFUN(Y,FY,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the matrix J and
%   the flag FLAG, the DJACFUN function must also set NEW_DATA. Otherwise, 
%   it should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead
%   to unnecessary copying).
%
%   The function DJACFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also KINSetOptions
%
%   NOTE: DJACFUN is specified through the property JacobianFn to KINSetOptions 
%   and is used only if the property LinearSolver was set to 'Dense'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2011/06/01 21:01:10 $
