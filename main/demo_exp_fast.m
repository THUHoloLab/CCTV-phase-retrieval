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
clear;clc;
close all;

% load experimental data
group_num = 1;
img_bg  = im2double(rgb2gray(imread(['../data/experiment/E',num2str(group_num),'/bg.bmp'])));
img_obj = im2double(rgb2gray(imread(['../data/experiment/E',num2str(group_num),'/obj.bmp'])));
load(['../data/experiment/E',num2str(group_num),'/params.mat'])

% normalization of the hologram
y = img_obj./mean(img_bg(:));

% select area of interest for reconstruction
figure
[temp,rect] = imcrop(img_obj);
if rem(size(temp,1),2) == 1
    rect(4) = rect(4) - 1;
end
if rem(size(temp,2),2) == 1
    rect(3) = rect(3) - 1;
end
close
y = imcrop(y,rect);
[n1,n2] = size(y);

% calculation of padding sizes to avoid circular boundary artifact
kernelsize  = params.dist*params.wavlen/params.pxsize/2;
nullpixels = ceil(kernelsize / params.pxsize);

% pre-calculate the transfer functions for diffraction modeling
H_f = fftshift(transfunc(n2+nullpixels*2,n1+nullpixels*2, params.dist,params.pxsize,params.wavlen,params.method)); % forward propagation
H_b = fftshift(transfunc(n2+nullpixels*2,n1+nullpixels*2,-params.dist,params.pxsize,params.wavlen,params.method)); % backward propagation

% forward model
Q  = @(x) ifft2(fft2(x).*H_f);      % forward propagation
QH = @(x) ifft2(fft2(x).*H_b);      % Hermitian of Q: backward propagation
C  = @(x) imgcrop(x,nullpixels);    % image cropping operation (to model the finite size of the sensor area)
CT = @(x) zeropad(x,nullpixels);    % transpose of C: zero-padding operation
A  = @(x) C(Q(x));                  % overall sampling operation
AH = @(x) QH(CT(x));                % Hermitian of A

%%
% =========================================================================
% Compressive phase retrieval algorithm
% =========================================================================
clear functions     % release memory (if using puma)

gpu = false;        % whether using GPU or not

% define the constraint
constraint = 'a';           % 'none': no constraint, 'a': absorption constraint only, 
                            % 's': support constraint only, 'as': absorption + support constraints
absorption = 1.1;           % define the upper bound for the modulus
support = zeros(n1+nullpixels*2,n2+nullpixels*2);   % define the support region
support(nullpixels+1:nullpixels+n1,nullpixels+1:nullpixels+n2) = 1;

% algorithm settings
x_est = AH(sqrt(y));    % initial guess
lam = 1e-2;             % regularization parameter
gam = 2;                % step size
n_iters    = 500;       % number of iterations (main loop)
n_subiters = 7;         % number of iterations (denoising)

% auxilary variables
z_est = x_est;
v_est = zeros(size(x_est,1),size(x_est,2),2);
w_est = zeros(size(x_est,1),size(x_est,2),2);

% initialize GPU
if gpu
    device  = gpuDevice();
    x_est   = gpuArray(x_est);
    y       = gpuArray(y);
    H_f     = gpuArray(H_f);
    H_b     = gpuArray(H_b);
    support = gpuArray(support);
    z_est   = gpuArray(z_est);
    v_est   = gpuArray(v_est);
    w_est   = gpuArray(w_est);
end

% define the fast projection operator
if strcmpi(constraint,'none')
    projf = @(x) x;
elseif strcmpi(constraint,'a')
    projf = @(x) min(abs(x),absorption).*exp(1i*angle(x));
elseif strcmpi(constraint,'s')
    projf = @(x) x.*support;
elseif strcmpi(constraint,'as')
    projf = @(x) min(abs(x),absorption).*exp(1i*angle(x)).*support;
else
    error("Invalid constraint. Should be 'A'(absorption), 'S'(support), 'AS'(both), or 'none'.")
end

% main loop
timer = tic;
for iter = 1:n_iters

    % print status
    fprintf('iter: %4d / %4d \n', iter, n_iters);
    
    % gradient update
    u = A(z_est);
    u = 1/2 * AH((abs(u) - sqrt(y)) .* exp(1i*angle(u)));
    u = z_est - gam * u;
    
    % proximal update
    v_est(:) = 0; w_est(:) = 0;
    for subiter = 1:n_subiters
        w_next = v_est + 1/8/gam*Df(projf(u-gam*DTf(v_est)));
        w_next = min(abs(w_next),lam).*exp(1i*angle(w_next));
        v_est = w_next + subiter/(subiter+3)*(w_next-w_est);
        w_est = w_next;
    end
    x_next = projf(u - gam*DTf(w_est));
    
    % Nesterov extrapolation
    z_est = x_next + (iter/(iter+3))*(x_next - x_est);
    x_est = x_next;
end

% wait for GPU
if gpu; wait(device); end
toc(timer)

% gather data from GPU
if gpu
    x_est   = gather(x_est);
    y       = gather(y);
    H_f     = gather(H_f);
    H_b     = gather(H_b);
    support = gather(support);
%     reset(device);
end

%%
% =========================================================================
% Display results
% =========================================================================
x_crop = x_est(nullpixels+1:nullpixels+n1,nullpixels+1:nullpixels+n2);
amp_est = abs(x_crop);
pha_est = puma_ho(angle(x_crop),1);

% visualize the reconstructed image
figure
set(gcf,'unit','normalized','position',[0.2,0.3,0.6,0.4])
subplot(1,2,1),imshow(amp_est,[]);colorbar
title('Retrieved amplitude','interpreter','latex','fontsize',14)
subplot(1,2,2),imshow(pha_est,[]);colorbar
title('Retrieved phase','interpreter','latex','fontsize',14)

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


function H = transfunc(nx, ny, dist, pxsize, wavlen, method)
% =========================================================================
% Calculate the transfer function of the free-space diffraction.
% -------------------------------------------------------------------------
% Input:    - nx, ny   : The image dimensions.
%           - dist     : Propagation distance.
%           - pxsize   : Pixel (sampling) size.
%           - wavlen   : Wavelength of the light.
%           - method   : Numerical method ('Fresnel' or 'Angular Spectrum').
% Output:   - H        : Transfer function.
% =========================================================================

% sampling in the spatial frequency domain
kx = pi/pxsize*(-1:2/nx:1-2/nx);
ky = pi/pxsize*(-1:2/ny:1-2/ny);
[KX,KY] = meshgrid(kx,ky);

k = 2*pi/wavlen;    % wave number

ind = (KX.^2 + KY.^2 >= k^2);  % remove evanescent orders
KX(ind) = 0; KY(ind) = 0;

if strcmp(method,'Fresnel')
    H = exp(1i*k*dist)*exp(-1i*dist*(KX.^2+KY.^2)/2/k);
elseif strcmp(method,'Angular Spectrum')
    H = exp(1i*dist*sqrt(k^2-KX.^2-KY.^2));
else
    errordlg('Wrong parameter for [method]: must be <Angular Spectrum> or <Fresnel>','Error');
end
end


function w = Df(x)
% =========================================================================
% Calculate the 2D gradient (finite difference) of an input image.
% -------------------------------------------------------------------------
% Input:    - x  : The input 2D image.
% Output:   - w  : The gradient (3D array).
% =========================================================================
w = cat(3,x(1:end,:) - x([2:end,end],:),x(:,1:end) - x(:,[2:end,end]));
end


function u = DTf(w)
% =========================================================================
% Calculate the transpose of the gradient operator.
% -------------------------------------------------------------------------
% Input:    - w  : 3D array.
% Output:   - x  : 2D array.
% =========================================================================
u1 = w(:,:,1) - w([end,1:end-1],:,1);
u1(1,:) = w(1,:,1);
u1(end,:) = -w(end-1,:,1);

u2 = w(:,:,2) - w(:,[end,1:end-1],2);
u2(:,1) = w(:,1,2);
u2(:,end) = -w(:,end-1,2);

u = u1 + u2;
end
