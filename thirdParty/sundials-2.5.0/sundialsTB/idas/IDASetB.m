function status = IDASetB(idxB, varargin)
%IDASetB changes optional input values during the integration.
%
%   Usage: IDASetB( IDXB, 'NAME1',VALUE1,'NAME2',VALUE2,... )
%
%   IDASetB can be used to change some of the optional inputs for
%   the backward problem identified by IDXB during the backward
%   integration, i.e., without need for a solver reinitialization.
%   The property names accepted by IDASet are a subset of those valid
%   for IDASetOptions. Any unspecified properties are left unchanged.
%   
%   IDASetB with no input arguments displays all property names.
%
%IDASetB properties
%(See also the IDAS User Guide)
%
%UserData - problem data passed unmodified to all user functions.
%  Set VALUE to be the new user data.  
%RelTol - Relative tolerance
%  Set VALUE to the new relative tolerance
%AbsTol - absolute tolerance
%  Set VALUE to be either the new scalar absolute tolerance or
%  a vector of absolute tolerances, one for each solution component.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/12/05 21:58:19 $

if (nargin == 0)
  fprintf('        UserData\n');
  fprintf('\n');
  fprintf('          RelTol\n');
  fprintf('          AbsTol\n');
  fprintf('\n');
  return;
end

KeyNames = {
    'UserData'
    'RelTol'
    'AbsTol'
           };

options = idm_options(KeyNames,varargin{:});

mode = 34;

status = idm(mode, idxB, options);
