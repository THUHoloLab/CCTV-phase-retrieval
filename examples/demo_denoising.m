% *************************************************************************
% * This code applies a fast gradient projection algorithm to
%   complex-valued image denoising problems.
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

arg = img + 0.1*randn(size(img));
mod = ones(size(img));
x = mod.*exp(1i*(arg*pi-pi/2));

figure
subplot(1,2,1),imshow(img*pi-pi/2,[-pi/2-0.5,pi/2+0.5]);colorbar;
title('Ground truth','interpreter','latex','fontsize',12)
subplot(1,2,2),imshow(angle(x),[-pi/2-0.5,pi/2+0.5]);colorbar;
title('Noisy observation','interpreter','latex','fontsize',12)
set(gcf,'unit','normalized','position',[0.2,0.3,0.6,0.4])

%% denoising
n_iters = 50;       % number of iteration for the denoising algorithm
gamma = 2e-1;       % regularization parameter

x_r = proxTV1i(c2r(x),gamma,n_iters);

figure,imshow(angle(r2c(x_r)),[-pi/2-0.5,pi/2+0.5]);colorbar
title(['Reconstruction after ',num2str(n_iters),' iterations'],'interpreter','latex','fontsize',12)
set(gcf,'unit','normalized','position',[0.3,0.3,0.4,0.4])