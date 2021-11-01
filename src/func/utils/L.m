function grad = L(x)

[~,n1,n2,n3] = size(x);
grad = zeros(2,n1,n2,n3,3);

u = squeeze(x(1,:,:,:));
v = squeeze(x(2,:,:,:));

grad(1,:,:,:,1) = u - circshift(u,[-1,0,0]);
grad(1,n1,:,:,1) = 0;
grad(1,:,:,:,2) = u - circshift(u,[0,-1,0]);
grad(1,:,n2,:,2) = 0;
grad(1,:,:,:,3) = u - circshift(u,[0,0,-1]);
grad(1,:,:,n3,3) = 0;

grad(2,:,:,:,1) = v - circshift(v,[-1,0,0]);
grad(2,n1,:,:,1) = 0;
grad(2,:,:,:,2) = v - circshift(v,[0,-1,0]);
grad(2,:,n2,:,2) = 0;
grad(2,:,:,:,3) = v - circshift(v,[0,0,-1]);
grad(2,:,:,n3,3) = 0;

end

