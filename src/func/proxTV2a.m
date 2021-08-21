function prox = proxTV2a(x,gamma,iter,alpha,varargin)

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
    if alpha > 0 && alpha < 1
        deno(1,:,:,:) = max(alpha,abs(grad_next(1,:,:,:)));
        deno(2,:,:,:) = max(alpha,abs(grad_next(2,:,:,:)));
        deno(3,:,:,:) = max(alpha,abs(grad_next(3,:,:,:)));
        deno(4,:,:,:) = max(1-alpha,abs(grad_next(4,:,:,:)));
        deno(5,:,:,:) = max(1-alpha,abs(grad_next(5,:,:,:)));
        deno(6,:,:,:) = max(1-alpha,abs(grad_next(6,:,:,:)));
        grad_next = grad_next./deno;
        grad_next(1:3,:,:,:) = grad_next(1:3,:,:,:)*alpha;
        grad_next(4:6,:,:,:) = grad_next(4:6,:,:,:)*(1-alpha);
    elseif alpha == 0
        deno(4,:,:,:) = max(1,abs(grad_next(4,:,:,:)));
        deno(5,:,:,:) = max(1,abs(grad_next(5,:,:,:)));
        deno(6,:,:,:) = max(1,abs(grad_next(6,:,:,:)));
        grad_next(1:3,:,:,:) = zeros(3,n1,n2,n3);
        grad_next(4:6,:,:,:) = grad_next(4:6,:,:,:)./deno(4:6,:,:,:);
    elseif alpha == 1
        deno(1,:,:,:) = max(1,abs(grad_next(1,:,:,:)));
        deno(2,:,:,:) = max(1,abs(grad_next(2,:,:,:)));
        deno(3,:,:,:) = max(1,abs(grad_next(3,:,:,:)));
        grad_next(4:6,:,:,:) = zeros(3,n1,n2,n3);
        grad_next(1:3,:,:,:) = grad_next(1:3,:,:,:)./deno(1:3,:,:,:);
    end
    t_next = (1+sqrt(1+4*t_prev^2))/2;
    temp = grad_next + (t_prev-1)/t_next*(grad_next-grad_prev);
    grad_prev = grad_next;
    t_prev = t_next;
end  

prox = proj(x - gamma*LT(grad_next));

end

