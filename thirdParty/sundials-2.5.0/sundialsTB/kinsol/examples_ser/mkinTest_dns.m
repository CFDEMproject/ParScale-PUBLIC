function mkinTest_dns
%mkinTest_dns - KINSOL example problem (serial, dense)
%   Simple test problem for the Dense linear solver in KINSOL
%   This example solves the system
%       y(1)^2 + y(2)^2 = 1
%       y(2) = y(1)^2
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/12/05 21:58:19 $

% Initialize problem
neq = 2;
fnormtol  = 1.0e-5;
scsteptol = 1.0e-4;
msbset = 1; % force exact Newton
options = KINSetOptions('FuncNormTol', fnormtol,...
                        'ScaledStepTol',scsteptol,...
                        'LinearSolver','Dense',....
                        'MaxNumSetups', msbset);
KINInit(@sysfn, neq, options);

% Solve problem
y0 = ones(neq,1);
scale = ones(neq,1);
strategy = 'LineSearch';
[status, y] = KINSol(y0, strategy, scale, scale);

% Evaluate system function at solution
[fy, flag] = sysfn(y);

% Print results
fprintf('Solution: %10.4e  %10.4e\n', y(1), y(2));
fprintf('Residual: %10.4e  %10.4e\n', fy(1), fy(2));

slv_stats = KINGetStats;
ls_stats = slv_stats.LSInfo;
slv_stats
ls_stats


% Free memory
KINFree;

return

% ===================================================================

function [fy, flag] = sysfn(y)

fy(1) = y(1)^2 + y(2)^2 - 1.0;
fy(2) = y(2) - y(1)^2;

flag = 0;

return


