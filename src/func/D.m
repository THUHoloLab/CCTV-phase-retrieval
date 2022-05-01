function w = D(x)

[n1,n2] = size(x);
w = zeros(n1,n2,2);

w(:,:,1) = x - circshift(x,[-1,0]);
w(n1,:,1) = 0;
w(:,:,2) = x - circshift(x,[0,-1]);
w(:,n2,2) = 0;

end

