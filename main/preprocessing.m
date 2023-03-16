% ========================================================================
% Introduction
% ========================================================================
% This code shows the preprocessing steps needed to run the phase retrieval
% algorithm using your own experimental data.
%
% Reference:
%   - Y. Gao and L. Cao, ¡°Iterative projection meets sparsity 
%     regularization: towards practical single-shot quantitative phase 
%     imaging with in-line holography,¡± Light: Advanced Manufacturing 4, 
%     1-17 (2023).
% 
% Author: Yunhui Gao (gyh21@mails.tsinghua.edu.cn)
% =========================================================================
%% generate data
clear;clc;
close all;

% load functions and data
addpath(genpath('../src'))
addpath(genpath('./utils'))
group_num = 1;
img_bg  = im2double(rgb2gray(imread(['../data/experiment/E',num2str(group_num),'/bg.bmp'])));   % reference image
img_obj = im2double(rgb2gray(imread(['../data/experiment/E',num2str(group_num),'/obj.bmp'])));  % hologram

y = img_obj./mean(img_bg(:));   % divide the hologram by a reference image

% select area of calculation
disp('Please select the area for diffraction calculation ...')
figure
[temp,rect_cal] = imcrop(y);
if rem(size(temp,1),2) == 1
    rect_cal(4) = rect_cal(4) - 1;
end
if rem(size(temp,2),2) == 1
    rect_cal(3) = rect_cal(3) - 1;
end
close
disp('Selected.')
y_crop = imcrop(y,rect_cal);

% select area of display
disp('Please select the area for display ...')
figure
[temp,rect_dis] = imcrop(y_crop);
if rem(size(temp,1),2) == 1
    rect_dis(4) = rect_dis(4) - 1;
end
if rem(size(temp,2),2) == 1
    rect_dis(3) = rect_dis(3) - 1;
end
close
disp('Selected.')

% physical parameters
params.pxsize = 5.86e-3;                % pixel size (mm)
params.wavlen = 0.660e-3;               % wavelength (mm)
params.method = 'Angular Spectrum';     % numerical method

%% determine the optimal sample-to-sensor distance
bias = 0.02;    figw = 0.50;    figh = 0.40;
figure,set(gcf,'unit','normalized','position',[(1-figw)/2,(1-figh)/2,figw,figh],'color','w')
[~, pos] = tight_subplot(1,2,[bias bias],[bias bias+0.04],[bias bias]);
for d = 5.4
    wavefront = propagate(sqrt(y_crop),-d, params.pxsize, params.wavlen, params.method);
    ax = subplot(1,2,1);
    ax.Position = pos{1};
    imshow(imcrop(abs(wavefront), rect_dis),[],'border','tight');
    title(['Amplitude (d = ',num2str(d),' mm)'])
    ax = subplot(1,2,2);
    ax.Position = pos{2};
    imshow(imcrop(angle(wavefront), rect_dis),[],'border','tight');
    title(['Phase (d = ',num2str(d),' mm)'])
    drawnow;
    pause(1)
end

%% set the distance
params.dist = 5.4;      % sample-to-sensor distance (mm)

%% save the calibrated data
save(['../data/experiment/E',num2str(group_num),'/params.mat'],'params')