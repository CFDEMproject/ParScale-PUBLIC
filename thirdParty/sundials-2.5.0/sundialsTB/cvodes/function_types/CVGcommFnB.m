%CVGcommFn - type for user provided communication function (BBDPre) for backward problems.
%
%   The function GCOMFUNB must be defined either as
%        FUNCTION FLAG = GCOMFUNB(T, Y, YB)
%   or as
%        FUNCTION [FLAG, NEW_DATA] = GCOMFUNB(T, Y, YB, DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. 
%
%   The function GCOMFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVGlocalFnB, CVodeSetOptions
%
%   NOTES:
%     GCOMFUNB is specified through the GcommFn property in CVodeSetOptions
%     and is used only if the property PrecModule is set to 'BBDPre'.
%
%     Each call to GCOMFUNB is preceded by a call to the RHS function
%     ODEFUNB with the same arguments T, Y, and YB. Thus GCOMFUNB can
%     omit any communication done by ODEFUNB if relevant to the evaluation
%     of G by GLOCFUNB. If all necessary communication was done by ODEFUNB,
%     GCOMFUNB need not be provided.     

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2011/06/01 20:44:05 $
