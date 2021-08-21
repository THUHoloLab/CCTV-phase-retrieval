function c = r2c(r)

% c = r(:,:,1) + 1i*r(:,:,2);

% [n,~] = size(r);
% c = r(1:n/2,:) + 1i*r(n/2+1:n,:);

c = squeeze(r(1,:,:)) + 1i*squeeze(r(2,:,:));
    
end