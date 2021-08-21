function [x,n_iters,J_vals,runtimes] = FISTA(F,dF,lambda,x_init,varargin)
% *************************************************************************
% * This function applies the fast iterative shrinkage / thresholding 
%   algorithm (FISTA) to solve linear inverse problem of the form: 
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
%   [1] A. Beck and M. Teboulle, "A Fast Iterative Shrinkage-Thresholding 
%       Algorithm for Linear Inverse Problems," SIAM Journal on Imaging 
%       Sciences 2, 183-202 (2009).
%
% *************************************************************************
% * Author : Yunhui Gao
% * Date   : 2021/04/20
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
%   - 'AT'       : function handle, default = A'
%                  The transpose of A. Required when A is a
%                  function handle.
% 
%   - 'max_iter' : int, default = 100
%                  Maximum number of iterations allowed.
%
%   - 'min_iter' : int, default = 0
%                  Minimum number of iterations allowed.
%
%   - 'initializer'    : must be one of {0,1,2,array}, default = 0
%                        Initializer for x:
%                        0 -> initialization at zero,
%                        1 -> random initialization,
%                        2 -> initialization with AT(y),
%                        array -> initialization provided by the user
%                                 (same size as x).
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
%   - 'monotone' : bool, default = false
%                  Enforces monotonic decrease in J.
%
%   - 'verbose'  : bool, default = false
%                  Visualizes the iterative process.
%
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
r = @norml1;
proxr = @proxl1;
max_iter = 100;
min_iter = 0;
init = 0;
stop_criterion = 0;
tol = 1e-2;
eta = 2;
Lip = 1;
monotone = false;
verbose = false;

%% parse input arguments
for i = 1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'penalty'
            r = varargin{i+1};
        case 'prox_op'
            proxr = varargin{i+1};
        case 'at'
            AT = varargin{i+1};
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
        case 'monotone'
            monotone = varargin{i+1};
        case 'verbose'
            verbose = varargin{i+1};
        otherwise
            error(['Invalid parameter: ',varargin{i}]);
    end
end

if (sum(stop_criterion == [0 1 2 3]) == 0)
	error('Unknwon stopping criterion (''stop_criterion'')');
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
function val = Q(x,y,L)
    val = F(y) + dotArr(x-y,dF(y)) + L/2*normArr(x-y)^2 + R(x);
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
J_vals(1) = J(x);
L = Lip;
iter = 1;
loop = true;
crit = Inf;
t = 1;
z = x;
timer = tic;
while loop
    % line search
    while true 
        p = proxR(z - 1/L*dF(z), 1/L);
        if J(p) <= Q(p,z,L) && Q(p,z,L) <= Q(z,z,L)
            break
        end
        L = L * eta;
    end
    
    t_next = (1+sqrt(1+4*t^2))/2;
    
	if monotone 
        if J(p) <= J(x)
            x_next = p; 
        end
        z = x_next + (t-1)/t_next*(x_next - x) ...
            + t/t_next*(p - x_next);
    else
        x_next = p;
        z = x_next + (t-1)/t_next*(x_next - x);
	end

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
            iter, 'FISTA', J_vals(iter+1), 1/L, crit, runtimes(iter));
    end
    
    
    iter = iter + 1;
    x = x_next;
    t = t_next;
   
    
end

% display result
if verbose
    fprintf('Algorithm terminated.\n')
end

runtimes(isnan(runtimes)) = [];
J_vals(isnan(J_vals)) = [];
n_iters = iter - 1;
    
end