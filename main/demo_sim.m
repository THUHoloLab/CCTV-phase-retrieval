% ========================================================================
% Introduction
% ========================================================================
% This code provides a simple demonstration of compressive phase retrieval
% via constrained complex total variation (CCTV) regularization.
%
% Reference:
%   - Y. Gao and L. Cao, "Iterative projection meets sparsity 
%     regularization: towards practical single-shot quantitative phase 
%     imaging with in-line holography," Light: Advanced Manufacturing 4, 
%     6 (2023).
% 
% Author: Yunhui Gao (gyh21@mails.tsinghua.edu.cn)
% =========================================================================
%%
% =========================================================================
% Data generation
% =========================================================================
clear;clc
close all

% load functions
addpath(genpath('./utils'))
addpath(genpath('../src'))

% load test image
n = 256;
img = imresize(im2double(imread('../data/simulation/cameraman.bmp')),[n,n]);

% sample
x = exp(1i*pi*img);

% physical parameters
params.pxsize = 5e-3;                   % pixel size (mm)
params.wavlen = 0.5e-3;                 % wavelength (mm)
params.method = 'Angular Spectrum';     % numerical method
params.dist   = 5;                      % imaging distance (mm)

% zero-pad the object to avoid convolution artifacts
kernelsize = params.dist*params.wavlen/params.pxsize/2; % diffraction kernel size
nullpixels = ceil(kernelsize / params.pxsize);          % number of padding pixels
x = zeropad(x,nullpixels);                              % zero-padded sample

% forward model
Q  = @(x) propagate(x,params.dist,params.pxsize,params.wavlen,params.method);   % forward propagation
QH = @(x) propagate(x,-params.dist,params.pxsize,params.wavlen,params.method);  % Hermitian of Q: backward propagation
C  = @(x) imgcrop(x,nullpixels);    % image cropping operation (to model the finite size of the sensor area)
CT = @(x) zeropad(x,nullpixels);    % transpose of C: zero-padding operation
A  = @(x) C(Q(x));                  % overall sampling operation
AH = @(x) QH(CT(x));                % Hermitian of A

% generate data
rng(0)              % random seed, for reproducibility
noisevar = 0.01;    % noise level
y = abs(A(x)).^2;
y = y.*(1 + noisevar*randn(size(y)));      % Gaussian noise

% display
figure
set(gcf,'unit','normalized','Position',[0.15,0.3,0.7,0.4])
subplot(1,3,1),imshow(abs(x),[]);colorbar;
title('Amplitude of the object','interpreter','latex','fontsize',12)
subplot(1,3,2),imshow(angle(x),[]);colorbar;
title('Phase of the object','interpreter','latex','fontsize',12)
subplot(1,3,3),imshow(y,[]);colorbar;
title('Intensity measurement','interpreter','latex','fontsize',12)

%%
% =========================================================================
% Compressive phase retrieval algorithm
% =========================================================================

clear functions     % release memory (if using puma)

% define the constraint
global constraint
constraint = 'as';  % 'none': no constraint, 'a': absorption constraint only, 
                    % 's': support constraint only, 'as': absorption + support constraints
global absorption   % define the upper bound for the modulus
absorption = 1;
global support      % define the support region
support = zeros(size(x));
support(nullpixels+1:nullpixels+n,nullpixels+1:nullpixels+n) = 1;

% region for computing the errors
region.x1 = nullpixels+1;
region.x2 = nullpixels+n;
region.y1 = nullpixels+1;
region.y2 = nullpixels+n;

% algorithm settings
x_init = AH(sqrt(y));   % initial guess
lam = 2e-3;             % regularization parameter
gam = 2;                % step size (see the paper for details)
n_iters    = 500;       % number of iterations (main loop)
n_subiters = 1;         % number of iterations (denoising)

% options
opts.verbose = true;                                % print status during the iterations
opts.errfunc = @(z) relative_error_2d(z,x,region);  % user-defined error metrics
opts.display = true;                                % display intermediate results during the iterations
opts.autosave = false;                              % save the intermediate results

myF     = @(x) F(x,y,A);                        % fidelity function 
mydF    = @(x) dF(x,y,A,AH);                    % gradient of the fidelity function
myR     = @(x) CCTV(x,lam);                     % regularization function
myproxR = @(x,gam) prox(x,gam,lam,n_subiters);  % proximal operator for the regularization function

% run the algorithm
[x_est,J_vals,E_vals,runtimes] = APG(x_init,myF,mydF,myR,myproxR,gam,n_iters,opts);     % accelerated proximal gradient algorithm

%%
% =========================================================================
% Display results
% =========================================================================

% crop image to match the size of the sensor
x_crop = x_est(nullpixels+1:nullpixels+n,nullpixels+1:nullpixels+n);

% visualize the reconstructed image
figure
set(gcf,'unit','normalized','position',[0.2,0.3,0.6,0.4])
subplot(1,2,1),imshow(abs(x_crop),[]);colorbar
title('Retrieved amplitude','interpreter','latex','fontsize',14)
subplot(1,2,2),imshow(angle(x_crop),[]);colorbar
title('Retrieved phase','interpreter','latex','fontsize',14)

% visualize the convergence curves
figure,plot(0:n_iters,E_vals,'linewidth',1.5);
xlabel('Iteration')
ylabel('Error metric')
set(gca,'fontsize',14)

figure,plot(0:n_iters,J_vals,'linewidth',1.5);
xlabel('Iteration')
ylabel('Objective value')
set(gca,'fontsize',14)

%%
% =========================================================================
% Auxiliary functions
% =========================================================================

function v = F(x,y,A)
% =========================================================================
% Data-fidelity function.
% -------------------------------------------------------------------------
% Input:    - x   : The complex-valued transmittance of the sample.
%           - y   : Intensity image.
%           - A   : The sampling operator.
% Output:   - v   : Value of the fidelity function.
% =========================================================================
v = 1/2 * norm2(abs(A(x)) - sqrt(y))^2;

function n = norm2(x)   % calculate the l2 vector norm
n = norm(x(:),2);
end

end


function g = dF(x,y,A,AH)
% =========================================================================
% Gradient of the data-fidelity function.
% -------------------------------------------------------------------------
% Input:    - x   : The complex-valued transmittance of the sample.
%           - y   : Intensity image.
%           - A   : The sampling operator.
%           - AH  : Hermitian of A.
% Output:   - g   : Wirtinger gradient.
% =========================================================================
u = A(x);
u = (abs(u) - sqrt(y)) .* exp(1i*angle(u));
g = 1/2 * AH(u);
end


function u = imgcrop(x,cropsize)
% =========================================================================
% Crop the central part of the image.
% -------------------------------------------------------------------------
% Input:    - x        : Original image.
%           - cropsize : Cropping pixel number along each dimension.
% Output:   - u        : Cropped image.
% =========================================================================
u = x(cropsize+1:end-cropsize,cropsize+1:end-cropsize);
end


function u = zeropad(x,padsize)
% =========================================================================
% Zero-pad the image.
% -------------------------------------------------------------------------
% Input:    - x        : Original image.
%           - padsize  : Padding pixel number along each dimension.
% Output:   - u        : Zero-padded image.
% =========================================================================
u = padarray(x,[padsize,padsize],0);
end