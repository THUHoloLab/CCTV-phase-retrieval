function r = r2c_hermitian(c)

r = zeros([2,size(c)]);
r(1,:,:) = c;
r(2,:,:) = -1i*c;


end

