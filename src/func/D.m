function w = D(x)
% =========================================================================
% Calculate the 2D gradient (finite difference) of an input image.
% -------------------------------------------------------------------------
% Input:    - x  : 2D image.
% Output:   - w  : The gradient (3D array).
% =========================================================================

[n1,n2] = size(x);
w = zeros(n1,n2,2);

w(:,:,1) = x - circshift(x,[-1,0]);
w(n1,:,1) = 0;
w(:,:,2) = x - circshift(x,[0,-1]);
w(:,n2,2) = 0;

end

