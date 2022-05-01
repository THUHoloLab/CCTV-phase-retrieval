%% generate data
clear;clc;
% close all;

% load functions and test image
addpath(genpath('./utils'))
addpath(genpath('../src'))
img1 = im2double(imread('../data/simulation/cameraman.bmp'));
img2 = im2double(imread('../data/simulation/peppers.bmp'));

N = 256;
img1 = imresize(img1,[N,N]);
img2 = imresize(img2,[N,N]);

% object
x = exp(1i*pi*img1);

% physical parameters
params.pxsize = 5e-3;       % pixel size (mm)
params.wavlen = 0.5e-3;     % wavelength (mm)
params.method = 'Angular Spectrum';     % numerical method
params.dist = 5;            % imaging distance (mm)

dist_crit = 2*max([size(x,1),size(x,2)])*params.pxsize^2/max(params.wavlen);
if dist_crit < params.dist
    error('Angular spectrum not applicable')
end

kernelsize  = params.dist*params.wavlen/params.pxsize/2;
nullpixels = ceil(kernelsize / params.pxsize);

x = zeropad(x,nullpixels);

% forward model
Q  = @(x) propagate(x,params.dist,params.pxsize,params.wavlen,params.method);
QH = @(x) propagate(x,-params.dist,params.pxsize,params.wavlen,params.method);
C  = @(x) imgcrop(x,nullpixels);
CT = @(x) zeropad(x,nullpixels);
A  = @(x) C(Q(x));
AH = @(x) QH(CT(x));

% generate data
rng(0)
noisevar = 0.01;
y = abs(A(x)).^2;
y = y.*(1 + noisevar*randn(size(y)));      % intensity (gaussian noise)

% display measurement
figure
subplot(1,2,1),imshow(angle(x),[]);colorbar;
title('Phase of the object','interpreter','latex','fontsize',12)
subplot(1,2,2),imshow(y,[]);colorbar;
title('Intensity measurement','interpreter','latex','fontsize',12)
set(gcf,'unit','normalized','position',[0.2,0.3,0.6,0.4])

%% run the algorithm
global constraint
constraint = 'as';  % 'none', 'a', 's', or 'as'
global absorption
absorption = 1;
global support
support = zeros(size(x));
support(nullpixels+1:nullpixels+N,nullpixels+1:nullpixels+N) = 1;

% region for computing the errors
region.x1 = nullpixels+1;
region.x2 = nullpixels+N;
region.y1 = nullpixels+1;
region.y2 = nullpixels+N;

% algorithm settings
x_init = AH(sqrt(y));   % initial guess
lam = 2e-3;             % regularization parameter
gamma = 2;              % step size

n_iters    = 500;       % number of iterations (main loop)
n_subiters = 1;         % number of iterations (denoising)

% options
opts.verbose = true;
opts.errfunc = @(z) relative_error_2d(z,x,region);

myF     = @(x) F(x,y,A);
mydF    = @(x) dF(x,y,A,AH);
myR     = @(x) CCTV(x,lam);
myproxR = @(x,gamma) prox(x,gamma,lam,n_subiters);

% run the algorithm
[x_est,J_vals,E_vals,runtimes] = fista(x_init,myF,mydF,myR,myproxR,gamma,n_iters,opts);

%%
x_crop = x_est(nullpixels+1:nullpixels+N,nullpixels+1:nullpixels+N);
figure,imshow(abs(x_crop),[]);colorbar
figure,imshow(angle(x_crop),[]);colorbar

figure,plot(0:n_iters,E_vals);

%%
% =========================================================================
% Auxiliary functions
% =========================================================================

function v = F(x,y,A)
% =========================================================================
% Data-fidelity function.
% -------------------------------------------------------------------------
% Input:    - x : The complex-valued image.
% Output:   - v : Value of the fidelity function.
% =========================================================================
v = 1/2 * norm2(abs(A(x)) - sqrt(y))^2;

function n = norm2(x)
n = norm(x(:),2);
end

end

function g = dF(x,y,A,AH)
% =========================================================================
% Wirtinger gradient of the data-fidelity function.
% -------------------------------------------------------------------------
% Input:    - x : The 3D complex-valued spatiotemporal datacube.
% Output:   - g : Wirtinger gradient.
% =========================================================================
g = NaN(size(x));
    u = A(x);
    u = (abs(u) - sqrt(y)) .* exp(1i*angle(u));
    g = 1/2 * AH(u);
end


function u = imgcrop(x,cropsize)
u = x(cropsize+1:end-cropsize,cropsize+1:end-cropsize);
end

function u = zeropad(x,padsize)
u = padarray(x,[padsize,padsize],0);
end