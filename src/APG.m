function [x,J_vals,E_vals,runtimes] = APG(x_init,F,dF,R,proxR,gamma,n_iters,opts)
% =========================================================================
% The accelerated proximal gradient (APG) algorithm aiming to solve the
% optimization problems in the following form:
%               min J(x) = F(x) + R(x),
%                x
% where F(x) is differentiable and R(x) is non-differentiable.
% -------------------------------------------------------------------------
% Input:    - x_init   : Initial estimate.
%           - F        : Function handle for F(x).
%           - dF       : Function handle for the gradient of F(x).
%           - R        : Function handle for R(x),
%           - proxR    : Function handle for the proximal operator of R(x).
%           - gamma    : The step size.
%           - n_iters  : Number of iterations.
%           - opts     : Optional settings.
% Output:   - x        : The estimated solution x.
%           - J_vals   : The values of J(x) during the iterations.
%           - E_vals   : The values of the error function during the iterations.
%           - runtimes : The runtimes stored during the iterations.
% =========================================================================

% initialization
x = x_init;
z = x;
J_vals = NaN(n_iters+1,1);
E_vals = NaN(n_iters+1,1);
runtimes = NaN(n_iters,1);

J_vals(1) = F(x)+R(x);
if isfield(opts,'errfunc') && isa(opts.errfunc,'function_handle')
    E_vals(1) = opts.errfunc(z);
end

if isfield(opts,'display') && opts.display
    fig = figure;
    set(fig,'unit','normalized','position',[0.2,0.2,0.6,0.5])
end

if isfield(opts,'autosave') && opts.autosave
    foldername = 'cache';
    while exist(foldername)
        foldername = [foldername,'_new'];
    end
    mkdir(foldername)
end

timer = tic;
for iter = 1:n_iters
    
    % gradient projection update
    x_next = proxR(z - gamma*dF(z),gamma);
    J_vals(iter+1) = F(x_next)+R(x_next);
    z = x_next + (iter/(iter+3))*(x_next - x);
    
    % record runtime
    runtimes(iter) = toc(timer);
    
    % calculate error metric
    if isfield(opts,'errfunc') && isa(opts.errfunc,'function_handle')
        E_vals(iter+1) = opts.errfunc(z);
    end
    
    % print status
    if isfield(opts,'verbose') && opts.verbose
        fprintf('iter: %4d | objective: %10.4e | stepsize: %2.2e | runtime: %5.1f s\n', ...
                iter, J_vals(iter+1), gamma, runtimes(iter));
    end

    if isfield(opts,'autosave') && opts.autosave
        save([foldername,'/iter_',num2str(iter),'.mat'],'x');
    end
    
    x = x_next;
    
    % display intermediate results
    if isfield(opts,'display') && opts.display
        subplot(1,2,1),imshow(abs(x(:,:,1)),[]);colorbar
        subplot(1,2,2),imshow(angle(x(:,:,1)),[]);colorbar
        drawnow;
    end
end

end

