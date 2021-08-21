function grad = L(x)

[~,n1,n2,n3] = size(x);
grad = zeros(6,n1,n2,n3);

u = squeeze(x(1,:,:,:));
v = squeeze(x(2,:,:,:));

grad(1,:,:,:) = u - circshift(u,[-1,0,0]);
grad(1,n1,:,:) = 0;
grad(2,:,:,:) = u - circshift(u,[0,-1,0]);
grad(2,:,n2,:) = 0;
grad(3,:,:,:) = u - circshift(u,[0,0,-1]);
grad(3,:,:,n3) = 0;

grad(4,:,:,:) = v - circshift(v,[-1,0,0]);
grad(4,n1,:,:) = 0;
grad(5,:,:,:) = v - circshift(v,[0,-1,0]);
grad(5,:,n2,:) = 0;
grad(6,:,:,:) = v - circshift(v,[0,0,-1]);
grad(6,:,:,n3) = 0;

end

