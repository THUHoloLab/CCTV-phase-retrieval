function r = c2r(c)

% r = zeros([size(c),2]);
% r(:,:,1) = real(c);
% r(:,:,2) = imag(c);

% r = [real(c);imag(c)];

r = zeros([2,size(c)]);
r(1,:,:) = real(c);
r(2,:,:) = imag(c);

end

