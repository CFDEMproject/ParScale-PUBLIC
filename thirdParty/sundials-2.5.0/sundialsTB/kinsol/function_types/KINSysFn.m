%KINSysFn - type for user provided system function
%
%   The function SYSFUN must be defined as 
%        FUNCTION [FY, FLAG] = SYSFUN(Y)
%   and must return a vector FY corresponding to f(y).
%   If a user data structure DATA was specified in KINInit, then
%   SYSFUN must be defined as
%        FUNCTION [FY, FLAG, NEW_DATA] = SYSFUN(Y,DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then, besides setting the vector FY,
%   the SYSFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function SYSFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also KINInit
%
%   NOTE: SYSFUN is specified through the KINInit function.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2011/06/01 21:01:10 $
