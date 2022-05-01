function erg = energy_ho(kappa,psi,base,p,cliques,disc_bar,th,quant)
%energy_ho   Energy from kappa labeling and psi phase measurements.
%   erg = energy_ho(kappa,psi,base,p,cliques,disc_bar,p,th,quant) returns the energy of kappa labeling given the 
%   psi measurements image, the base ROI image (having ones in the region of interest (psi) and a passe-partout
%   made of zeros), the exponent p, the cliques matrix (each row indicating a displacement vector corresponding
%   to each clique), the disc_bar (complement to one of the quality maps), a threshold th defining a region for
%   which the potential (before a possible quantization) is quadratic, and quant which is a flag defining whether
%   the potential is or is not quantized.
%   (see J. Bioucas-Dias and G. Valadão, "Phase Unwrapping via Graph Cuts"
%   submitted to IEEE Transactions Image Processing, October, 2005).
%   SITE: www.lx.it.pt/~bioucas/ 

[m,n] = size(psi);
[cliquesm,cliquesn] = size(cliques); % Size of input cliques
maxdesl = max(max(abs(cliques))); % This is the maximum clique length used
% Here we put a passe-partout (constant length = maxdesl+1) in the images kappa and psi
base_kappa    = zeros(2*maxdesl+2+m,2*maxdesl+2+n); base_kappa(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1) = kappa;
psi_base      = zeros(2*maxdesl+2+m,2*maxdesl+2+n); psi_base(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1) = psi;
z = size(disc_bar,3);
base_disc_bar  = repmat(zeros(2*maxdesl+2+m,2*maxdesl+2+n),[1 1 z]); base_disc_bar(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,:) = disc_bar;

for t = 1:cliquesm
    % The allowed start and end pixels of the "interpixel" directed edge
    base_start(:,:,t) = circshift(base,[-cliques(t,1),-cliques(t,2)]).*base;
    base_end(:,:,t) = circshift(base,[cliques(t,1),cliques(t,2)]).*base;
    
    % By convention the difference images have the same size as the
    % original ones; the difference information is retrieved in the
    % pixel of the image that is subtracted (end of the diff vector)
    auxili = circshift(base_kappa,[cliques(t,1),cliques(t,2)]);
    t_dkappa(:,:,t) = (base_kappa-auxili);
    auxili2 = circshift(psi_base,[cliques(t,1),cliques(t,2)]);
    dpsi = auxili2 - psi_base;
    % Beyond base, we must multiply by
    % circshift(base,[cliques(t,1),cliques(t,2)]) in order to
    % account for frontier pixels that can't have links outside ROI
    a(:,:,t) = (2*pi*t_dkappa(:,:,t)-dpsi).*base.*circshift(base,[cliques(t,1),cliques(t,2)])...
               .*base_disc_bar(:,:,t);
end

erg = sum(sum(sum((clique_energy_ho(a,p,th,quant)))));

