% *************************************************************************
% * This code applies the proximal gradient method to phase retrieval.
% *************************************************************************
% * Author : Yunhui Gao
% * Date   : 2021/08/20
% *************************************************************************

%% generate data
clear;clc;
close all;

% load functions and test image
addpath(genpath('../src'))
addpath(genpath('./utils'))
img = im2double(imresize(imread('../data/cameraman.tif'),[256,256]));

% object
x = img;                  % pure amplitude object
x = exp(1i*(img)*pi*1/2);       % pure phase object

% physical parameters
params.pxsize = 5e-3;       % pixel size (mm)
params.wavlen = 0.5e-3;     % wavelength (mm)
params.method = 'Angular Spectrum';     % numerical method
params.dist = 10;           % imaging distance (mm)

% forward model
u = propagate(x,params.dist,params.pxsize,params.wavlen,params.method);
y = abs(u).^2.*(1 + 0.1*randn(size(u)));      % intensity

% display measurement
figure
subplot(1,2,1),imshow(angle(x),[]);colorbar;
title('Phase of the object','interpreter','latex','fontsize',12)
subplot(1,2,2),imshow(y,[]);colorbar;
title('Intensity measurement','interpreter','latex','fontsize',12)
set(gcf,'unit','normalized','position',[0.2,0.3,0.6,0.4])

%% define function handles
myF = @(x) F(x,y,params);       % the fidelity function F
mydF = @(x) dF(x,y,params);     % gradient of the fidelity function dF

%% initialization
x_init = propagate(sqrt(y),params.dist,params.pxsize,params.wavlen,params.method);
x_init = c2r(x_init);
x_init = proj(x_init);

%% run the algorithm
rng(0)              % random seed, for reproducibility
alpha = 0.5;
n_iters = 100;      % number of iterations for FISTA
n_subiters = 10;    % number of iterations to solve the denoising subproblem
lambda = 1e-2;

penalty = @(x) normTV1i(x,@indicator);    % isotropic TV norm as the penalty function
prox_op = @(x,gamma) proxTV1i(x,gamma,n_subiters,@proj);       % proximity operator

[x_r,~,~,~] = FISTA(myF,mydF,lambda,x_init,...        % run FISTA
'prox_op',      prox_op,...
'penalty',      penalty,...
'eta',          2,...
'Lip',          1,... 
'max_iter',     n_iters,...
'min_iter',     n_iters,...
'verbose',      true);

%% display results
x_c = r2c(x_r);
[phase_unwrap,~,~] = puma_ho(angle(x_c),1);     % phase unwrapping

figure
subplot(2,2,1),imshow(real(x_c),[]);colorbar
title('$\mathrm{Re}(\mathbf{x})$','interpreter','latex','fontsize',14)
subplot(2,2,2),imshow(imag(x_c),[]);colorbar
title('$\mathrm{Im}(\mathbf{x})$','interpreter','latex','fontsize',14)
subplot(2,2,3),imshow(abs(x_c),[]);colorbar
title('$\left| \mathbf{x} \right|$','interpreter','latex','fontsize',14)
subplot(2,2,4),imshow(phase_unwrap,[]);colorbar
title('$\mathrm{arg}(\mathbf{x})$','interpreter','latex','fontsize',14)
