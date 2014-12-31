%CVPrecSolveFnB - type for user provided preconditioner solve function for backward problems.
%
%   The user-supplied preconditioner solve function PSOLFN
%   is to solve a linear system P z = r in which the matrix P is
%   one of the preconditioner matrices P1 or P2, depending on the
%   type of preconditioning chosen.
%
%   The function PSOLFUNB must be defined either as
%        FUNCTION [ZB, FLAG] = PSOLFUNB(T,Y,YB,FYB,RB)
%   or as
%        FUNCTION [ZB, FLAG, NEW_DATA] = PSOLFUNB(T,Y,YB,FYB,RB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. In either case, it must return the vector ZB and the
%   flag FLAG.
%
%   See also CVPrecSetupFnB, CVodeSetOptions
%
%   NOTE: PSOLFUNB is specified through the property PrecSolveFn to 
%   CVodeSetOptions and is used only if the property LinearSolver was
%   set to 'GMRES', 'BiCGStab', or 'TFQMR' and if the property PrecType
%   is not 'None'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2011/06/01 20:44:05 $
