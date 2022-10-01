function [ w_o ] = propagate( w_i, dist, pxsize, wavlen, method )
% -----------------------------------------------------------------
% This function numerically simulates the free-space propagation of a 
% complex wavefield.
% -----------------------------------------------------------------
%   INPUT   w_i    : Input complex wavefield
%           dist   : Propagation distance
%           pxsize : Pixel size
%           wavlen : Wavelength
%           mtehod : Numerical method ('Fresnel' or 'Angular Spectrum')
%   OUTPUT  w_o    : Output wavefield after propagation
% -----------------------------------------------------------------                                                                                                                   function [ w_o ] = propagate( w_i, dist, pxsize, wavlen, method )

[N,M] = size(w_i);    % size of the wavefront

kx = pi/pxsize*(-1:2/M:1-2/M);
ky = pi/pxsize*(-1:2/N:1-2/N);
[KX,KY] = meshgrid(kx,ky);

k = 2*pi/wavlen;   % wave number

inputFT = fftshift(fft2(w_i));

if strcmp(method,'Fresnel')
    H = exp(1i*k*dist)*exp(-1i*dist*(KX.^2+KY.^2)/2/k);
elseif strcmp(method,'Angular Spectrum')
    % remove evanescent orders
    KX_m = KX;
    KY_m = KY;
    ind = (KX.^2+KY.^2 >= k^2);
    KX_m(ind) = 0;
    KY_m(ind) = 0;
    % transfer function
    H = exp(1i*dist*sqrt(k^2-KX_m.^2-KY_m.^2));
else
    errordlg('Wrong parameter for [method]: must be <Angular Spectrum> or <Fresnel>','Error');
end

w_o = ifft2(fftshift(inputFT.*H));

end

