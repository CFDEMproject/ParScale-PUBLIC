function status = CVodeSensToggleOff()
% CVodeSensToggleOff deactivates sensitivity calculations.
%   It does NOT deallocate sensitivity-related memory so that 
%   sensitivity computations can be later toggled ON (through
%   CVodeSensReInit).
%
%   Usage: CVodeSensToggleOff
%
%   See also: CVodeSensInit, CVodeSensReInit

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/12/05 21:58:18 $

mode = 18;
status = cvm(mode);
