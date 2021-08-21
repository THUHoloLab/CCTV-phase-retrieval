function prox = proxTV1a(x,gamma,iter,varargin)

if isempty(varargin)
    proj = @(x) x;
else
    proj = varargin{1};
end

t_prev = 1;

[~,n1,n2,n3] = size(x);
grad_next = zeros(6,n1,n2,n3);
grad_prev = zeros(6,n1,n2,n3);
temp = zeros(6,n1,n2,n3);

for i = 1:iter
    grad_next = temp + 1/8/gamma*L(proj(x - gamma*LT(temp)));
    deno = zeros(6,n1,n2,n3);
    deno(1,:,:,:) = max(1,sqrt(grad_next(1,:,:,:).^2 + grad_next(4,:,:,:).^2));
    deno(2,:,:,:) = max(1,sqrt(grad_next(2,:,:,:).^2 + grad_next(5,:,:,:).^2));
    deno(3,:,:,:) = max(1,sqrt(grad_next(3,:,:,:).^2 + grad_next(6,:,:,:).^2));
    deno(4,:,:,:) = deno(1,:,:,:);
    deno(5,:,:,:) = deno(2,:,:,:);
    deno(6,:,:,:) = deno(3,:,:,:);
    grad_next = grad_next./deno;
    t_next = (1+sqrt(1+4*t_prev^2))/2;
    temp = grad_next + (t_prev-1)/t_next*(grad_next-grad_prev);
    grad_prev = grad_next;
    t_prev = t_next;
end  

prox = proj(x - gamma*LT(grad_next));

end

