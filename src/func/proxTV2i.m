function prox = proxTV2i(x,gamma,iter,alpha,varargin)

if isempty(varargin)
    proj = @(x) x;
else
    proj = varargin{1};
end

t_prev = 1;

[~,n1,n2,n3] = size(x);
grad_next = zeros(2,n1,n2,n3,3);
grad_prev = zeros(2,n1,n2,n3,3);
temp = zeros(2,n1,n2,n3,3);

for i = 1:iter
    grad_next = temp + 1/8/gamma*L(proj(x - gamma*LT(temp)));
    deno = zeros(2,n1,n2,n3,3);
    if alpha > 0 && alpha < 1
        deno(1,:,:,:,1) = max(alpha,sqrt(grad_next(1,:,:,:,1).^2 + grad_next(1,:,:,:,1).^2 + grad_next(1,:,:,:,3).^2));
        deno(1,:,:,:,2) = deno(1,:,:,:,1);
        deno(1,:,:,:,3) = deno(1,:,:,:,1);
        deno(2,:,:,:,1) = max(1-alpha,sqrt(grad_next(2,:,:,:,1).^2 + grad_next(2,:,:,:,2).^2 + grad_next(2,:,:,:,3).^2));
        deno(2,:,:,:,2) = deno(2,:,:,:,1);
        deno(2,:,:,:,3) = deno(2,:,:,:,1);
        grad_next = grad_next./deno;
        grad_next(1,:,:,:,:) = grad_next(1,:,:,:,:)*alpha;
        grad_next(2,:,:,:,:) = grad_next(2,:,:,:,:)*(1-alpha);
    elseif alpha == 0
        deno(2,:,:,:,1) = max(1,sqrt(grad_next(2,:,:,:,1).^2 + grad_next(2,:,:,:,2).^2 + grad_next(2,:,:,:,3).^2));
        deno(2,:,:,:,2) = deno(2,:,:,:,1);
        deno(2,:,:,:,3) = deno(2,:,:,:,1);
        grad_next(1,:,:,:,:) = zeros(n1,n2,n3,3);
        grad_next(2,:,:,:,:) = grad_next(2,:,:,:,:)./deno(2,:,:,:,:);
    elseif alpha == 1
        deno(1,:,:,:,1) = max(1,sqrt(grad_next(1,:,:,:,1).^2 + grad_next(1,:,:,:,2).^2 + grad_next(1,:,:,:,3).^2));
        deno(1,:,:,:,2) = deno(1,:,:,:,1);
        deno(1,:,:,:,3) = deno(1,:,:,:,1);
        grad_next(2,:,:,:,:) = zeros(n1,n2,n3,3);
        grad_next(1,:,:,:,:) = grad_next(1,:,:,:,:)./deno(1,:,:,:,:);
    end
    t_next = (1+sqrt(1+4*t_prev^2))/2;
    temp = grad_next + (t_prev-1)/t_next*(grad_next-grad_prev);
    grad_prev = grad_next;
    t_prev = t_next;
end

prox = proj(x - gamma*LT(grad_next));

end

