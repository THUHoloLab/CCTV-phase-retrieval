function [x,n_iters,J_vals,errs,runtimes] = TwIST(F,dF,lambda,x_init,varargin)
% *************************************************************************
% * This function applies the two-Step iterative shrinkage / thresholding 
%   algorithm (TwIST) to solve linear inverse problem of the form: 
%
%           min {J(x) = F(x) + R(x)},
%            x
%   
%   where   F(x) = 0.5*|| A x - y ||_2^2    is the data-fidelity term,
%   and     R(x) = lambda * r(x)            is the regularization term.
%
%   See the Readme.md file for details.
% 
% * References:
%   [1] J. M. Bioucas-Dias and M. A. T. Figueiredo, "A New TwIST: Two-Step
%       Iterative Shrinkage/Thresholding Algorithms for Image Restoration,"
%       IEEE Transactions on Image Processing 16, 2992-3004 (2007).
%
% *************************************************************************
% * This implementation is based on the code by Jos¨¦ Bioucas-Dias and 
%   M¨¢rio Figueiredo (see http://www.lx.it.pt/~bioucas/TwIST/TwIST.htm).
%   -----------------------------------------------------------------------
% * Copyright (2010): Jos¨¦ Bioucas-Dias and M¨¢rio Figueiredo
%   TwIST is distributed under the terms of the GNU General Public
%   License 2.0.
% 
% * Permission to use, copy, modify, and distribute this software for any 
%   purpose without fee is hereby granted, provided that this entire
%   notice is included in all copies of any software which is or includes
%   a copy or modification of this software and in all copies of the
%   supporting documentation for such software.
%   This software is being provided "as is", without any express or
%   implied warranty.  In particular, the authors do not make any
%   representation or warranty of any kind concerning the merchantability
%   of this software or its fitness for any particular purpose."
% =========================================================================
% * Modified:   Yunhui Gao.
% * Date:       2021/04/20.
% *************************************************************************
%
%   ===== Required inputs =================================================
%
%	- y  : 1D / 2D / 3D array
%          Observation.
%
%	- A  : 2D array or function handle
%          The forward linear operator of the system. A can be a 2D array 
%          when both x and y are 1D arrays, otherwise A should be a 
%          function handle.
%
% 	- lambda : float
%              Regularization paramter.
%
%   ===== Optional inputs =================================================
%
%   - 'penalty'  : function handle, default = @norml1 
%                  The penalty function r(x), should accept one input x.
%                  The default choice is the l1 norm r(x) = || x ||_1.
%
%   - 'prox_op'  : function handle, default = @proxl1 ( soft thresholding )
%                  Proximity operator for r(x), should accept two inputs 
%                  x and gamma. The default choice is the soft-
%                  thresholding function.
%
%   - 'max_iter' : int, default = 100
%                  Maximum number of iterations allowed.
%
%   - 'min_iter' : int, default = 0
%                  Minimum number of iterations allowed.
%        
%   - 'stop_criterion' : must be one of {0,1,2,3}, default = 0
%                        Type of stopping criterion to use:
%                        0 -> stop when the relative change in the
%                             objective function falls below 'tol',
%                        1 -> stop when the relative norm of the  
%                             difference between two consecutive  
%                             estimates falls below 'tol',
%                        2 -> stop when the objective function 
%                             becomes equal or less than 'tol',
%                        3 -> stop only when the number of iteration
%                             reaches 'max_iter'.
%
%   - 'tol' : float, default = 0.01
%             Tolerance, used for stopping criterion.
%
%   - 'eta' : float, default = 2
%             A positive (>1) scalar used for backtracking.
%         
%   - 'lip' : int, default = 1
%             Esitimate of the Lipschitz constant, which helps determine
%             the stepsize.
%
%   - 'eig_min' : float, default = 1e-4
%                 Estimated minimum eigenvalue of A'*A. 
%                 If minimum eigenvalue of A'*A == 0, or unknwon,  
%                 set eig_min to a value much smaller than 1.
%
%                 Rule of Thumb: 
%                 -> eig_min = 1e-4 for severyly ill-conditioned problems
%                 -> eig_min = 1e-2 for mildly  ill-conditioned problems
%                 -> eig_min = 1    for A unitary direct operators
%
%   - 'alpha'  : float, default = calculated as a function of eig_min
%                Parameter alpha of TwIST (see Eq. (22) of ref. [1]).
%
%   - 'beta'   : float, default = calculated as a function of eig_min
%                Parameter beta of TwIST (see Eq. (23) of ref. [1]).
%
%   - 'monotone' : bool, default = true
%                  Enforces monotonic decrease in J.
%
%   - 'verbose'  : bool, default = false
%                  Visualizes the iterative process.
%
%   - 'ground_truth'  : 2D array, default: None
%                       The ground truth of x.
%
%   ===== Outputs =========================================================
%
%   - x  : 1D / 2D / 3D array
%          The solution.
%
%   - n_iters : int
%               Number of iterations.
%
%   - J_vals   : 1D array
%                The values of J(x) for each iteration.
%
%   - runtimes : 1D array
%                Runtime of the algorithm for each iteration.
%
% *************************************************************************

%% settings
% add path
addpath(genpath('func'));   % path for core functions
addpath(genpath('utils'));  % path for helper functions

% assign default values
r = @norml1;        % penalty function r(x)
proxr = @proxl1;    % proximity operator for the penalty function
max_iter = 100;
min_iter = 0;
stop_criterion = 0;
tol = 1e-2;
eta = 2;
Lip = 1;
eig_min = 1e-4;
eig_max = 1;
alpha = 0;
beta  = 0;
monotone = true;
verbose = false;
x_gt = NaN;

%% parse input arguments
for i = 1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'penalty'
            r = varargin{i+1};
        case 'prox_op'
            proxr = varargin{i+1};
        case 'max_iter'
            max_iter = varargin{i+1};
        case 'min_iter'
            min_iter = varargin{i+1};
        case 'stop_criterion'
            stop_criterion = varargin{i+1};
        case 'tol'
            tol = varargin{i+1};
        case 'eta'
            eta = varargin{i+1};
        case 'lip'
            Lip = varargin{i+1};
        case 'eig_min'
            eig_min = varargin{i+1};
        case 'alpha'
            alpha = varargin{i+1};
        case 'beta'
            beta = varargin{i+1};
        case 'monotone'
            monotone = varargin{i+1};
        case 'verbose'
            verbose = varargin{i+1};
        case 'ground_truth'
            x_gt = varargin{i+1};
        otherwise
            error(['Invalid parameter: ',varargin{i}]);
    end
end

if (sum(stop_criterion == [0 1 2 3]) == 0)
	error('Unknwon stopping criterion (''stop_criterion'')');
end

% set TwIST parameters 
rho0 = (1-eig_min/eig_max)/(1+eig_min/eig_max);
if alpha == 0 
    alpha = 2/(1+sqrt(1-rho0^2));
end
if beta == 0 
    beta  = alpha*2/(eig_min+eig_max);
end

%% 
% =========================================================================
%                         auxilary functions
% =========================================================================
% calculate the regularization term R(x)
function val = R(x)
    val = lambda*r(x);
end

% calculate the objective function J(x)
function val = J(x)
    val = F(x) + R(x);
end

% calculate the quadratic approximation function Q(x,y,L)
function val = Q(x,z,L)
    val = F(z) + dotArr(x-z,dF(z)) + L/2*normArr(x-z)^2 + R(x);
end

% calculate the proximity operator for function R(x)
function val = proxR(x,gamma)
    val = proxr(x,gamma*lambda);
end

%% 
% =========================================================================
%                               main loop
% =========================================================================
% initialization
x = x_init;

runtimes = NaN(max_iter,1);
J_vals = NaN(max_iter,1);
J_vals(1) = J(x);   % initial value for the objective function

L = Lip;        % estimate for the Lipschitz constant
iter = 1;       % counter for the iterations
crit = Inf;     % criterion for stopping
loop = true;    % boolean flag for the main loop
twist = false;  % boolean flag for the TwIST update

timer = tic;
while loop      % main loop
    
    while true  % inner loop (a single update)
        
        p = proxR(x - 1/L*dF(x), 1/L);      % IST iteration
        
        if twist                            
            x_est = (1-alpha)*x_prev + (alpha-beta)*x + beta*p;     % TwIST iteration
            
            if J(x_est) > J(x) && monotone  % enforce monotonicity
                twist = false;              % monotonicity violated, try IST instead
                
            else                            % inner loop finished, a TwIST update
                x_next = x_est;
                str = 'TwIST';      
                break;
            end
            
        else
            if Q(p,x,L) > Q(x,x,L) ...  % z = p should be the minimizer of Q(z,x,L) as a function of z.
                    ...                 % If not satisfied, may due to the inexact calculation of proxR.
                    || J(p) > Q(p,x,L)  % Assuming that p is exact, the Lipschitz estimate is not accurate.
                L = L * eta;            % increase the estimate for the Lip constant (decrease the stepsize)
                
            else                        % inner loop finished, an IST update
                twist = true;
                x_next = p;
                str = '  IST';
                break;
            end
        end
    end
    
    % record information
    J_vals(iter+1) = J(x_next);
    runtimes(iter) = toc(timer);
    
    % check the stopping criterion
    switch stop_criterion
        case 0
            crit = abs(J_vals(iter+1)-J_vals(iter))/J_vals(iter);
        case 1
            crit = normArr(x_next-x)/normArr(x);
        case 2
            crit = J_vals(iter+1);
    end
    if iter >= min_iter
        if iter >= max_iter || crit < tol
            loop = false;
        end
    end
    
    % display progress
    if verbose
        fprintf(['iter: %5d | type: %s | cost: %10.4e | stepsize: %2.2e | ', ...
            'criterion: %10.4e | runtime: %5.1f s\n'], ...
            iter, str, J_vals(iter+1), 1/L, crit, runtimes(iter));
    end
    
    % update for next iteration
    iter = iter + 1;
    x_prev = x;
    x = x_next;
    
end

% display result
if verbose
    fprintf('Algorithm terminated.\n')
end

runtimes(isnan(runtimes)) = [];
J_vals(isnan(J_vals)) = [];
n_iters = iter - 1;

end