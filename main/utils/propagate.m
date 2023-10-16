function [ w_o ] = propagate( w_i, dist, pxsize, wavlen, method )
% =========================================================================
% Calculate the free-space propagation of a complex wavefield.
% -------------------------------------------------------------------------
% Input:    - w_i    : Input complex wavefield
%           - dist   : Propagation distance
%           - pxsize : Pixel size
%           - wavlen : Wavelength
%           - mtehod : Numerical method ('Fresnel' or 'Angular Spectrum')
% Output:   - w_o    : Output wavefield after propagation
% =========================================================================

[nx,ny] = size(w_i);      % size of the wavefront

% sampling in the spatial frequency domain
kx = pi/pxsize*(-1:2/ny:1-2/ny);
ky = pi/pxsize*(-1:2/nx:1-2/nx);
[KX,KY] = meshgrid(kx,ky);

k = 2*pi/wavlen;        % wave number

% remove evanescent orders
ind = (KX.^2+KY.^2 >= k^2);
KX(ind) = 0; KY(ind) = 0;

inputFT = fftshift(fft2(w_i));

if strcmp(method,'Fresnel')
    H = exp(1i*k*dist)*exp(-1i*dist*(KX.^2+KY.^2)/2/k);
elseif strcmp(method,'Angular Spectrum')
    H = exp(1i*dist*sqrt(k^2-KX.^2-KY.^2));
else
    errordlg('Wrong parameter for [method]: must be <Angular Spectrum> or <Fresnel>','Error');
end

w_o = ifft2(fftshift(inputFT.*H));

end

