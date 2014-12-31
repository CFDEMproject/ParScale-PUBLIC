function status = IDASensReInit(yyS0,ypS0,options)
%IDASensReInit reinitializes IDAS's FSA-related memory
%   assuming it has already been allocated in prior calls to IDAInit
%   and IDASensInit.
%   The number of sensitivities Ns is assumed to be unchanged since the 
%   previous call to IDASensInit.
%
%   Usage: IDASensReInit ( YYS0, YPS0 [, OPTIONS ] ) 
%
%   YYS0, YPS0   Initial conditions for sensitivity variables.
%            YYS0 and YPS0 must be matrices with N rows and Ns columns, where N is 
%            the problem dimension and Ns the number of sensitivity systems.
%   OPTIONS  is an (optional) set of FSA options, created with
%            the IDASetFSAOptions function. 
%
%   See also: IDASensSetOptions, IDAReInit, IDASensInit

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/12/05 21:58:19 $

mode = 13;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

status = idm(mode, yyS0, ypS0, options);
