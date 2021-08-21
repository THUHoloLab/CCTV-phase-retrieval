function r = r2c_transpose(c)

% r = zeros([size(c),2]);
% r(:,:,1) = c;
% r(:,:,2) = 1i*c;

% r = [c;1i*c];

r = zeros([2,size(c)]);
r(1,:,:) = c;
r(2,:,:) = 1i*c;

end

