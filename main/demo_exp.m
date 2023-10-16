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

% load functions
addpath(genpath('./utils'))
addpath(genpath('../src'))

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

% forward model
Q  = @(x) propagate(x, params.dist,params.pxsize,params.wavlen,params.method);      % forward propagation
QH = @(x) propagate(x,-params.dist,params.pxsize,params.wavlen,params.method);      % Hermitian of Q: backward propagation
C  = @(x) imgcrop(x,nullpixels);                % image cropping operation (to model the finite size of the sensor area)
CT = @(x) zeropad(x,nullpixels);                % transpose of C: zero-padding operation
A  = @(x) C(Q(x));                              % overall sampling operation
AH = @(x) QH(CT(x));                            % Hermitian of A

%%
% =========================================================================
% Compressive phase retrieval algorithm
% =========================================================================
clear functions     % release memory (if using puma)

% define the constraint
global constraint
constraint = 'a';   % 'none': no constraint, 'a': absorption constraint only, 
                    % 's': support constraint only, 'as': absorption + support constraints
global absorption   % define the upper bound for the modulus
absorption = 1.1;
global support      % define the support region
support = zeros(n1+nullpixels*2,n2+nullpixels*2);
support(nullpixels+1:nullpixels+n1,nullpixels+1:nullpixels+n2) = 1;

% region for computing the errors
region.x1 = nullpixels+1;
region.x2 = nullpixels+n1;
region.y1 = nullpixels+1;
region.y2 = nullpixels+n2;

% algorithm settings
x_init = AH(sqrt(y));   % initial guess
lam = 1e-2;             % regularization parameter
gam = 2;                % step size
n_iters    = 500;       % number of iterations (main loop)
n_subiters = 7;         % number of iterations (denoising)

% options
opts.verbose = true;
opts.errfunc = nan;
opts.display = true;
opts.autosave = false;

% building block
myF     = @(x) F(x,y,A);                            % fidelity function 
mydF    = @(x) dF(x,y,A,AH);                        % gradient of the fidelity function
myR     = @(x) CCTV(x,lam);                         % regularization function
myproxR = @(x,gamma) prox(x,gamma,lam,n_subiters);  % proximal operator for the regularization function

% run the algorithm
[x_est,J_vals,E_vals,runtimes] = APG(x_init,myF,mydF,myR,myproxR,gam,n_iters,opts);   % accelerated proximal gradient algorithm


%%
% =========================================================================
% Display results
% =========================================================================
x_crop = x_est(nullpixels+1:nullpixels+n1,nullpixels+1:nullpixels+n2);
amp_est = abs(x_crop);
pha_est = puma_ho(angle(x_crop),1);

% visualize the reconstructed image
figure
subplot(1,2,1),imshow(amp_est,[]);colorbar
title('Retrieved amplitude','interpreter','latex','fontsize',14)
subplot(1,2,2),imshow(pha_est,[]);colorbar
title('Retrieved phase','interpreter','latex','fontsize',14)
set(gcf,'unit','normalized','position',[0.2,0.3,0.6,0.4])

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