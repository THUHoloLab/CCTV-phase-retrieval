function norm = normTV2a(x,alpha,varargin)

if isempty(varargin)
    indicator = @(x) 0;
else
    indicator = varargin{1};
end

grad = L(x);
norm = alpha*sum(sqrt(sum(grad(1,:,:,:,:).^2,5)),[2,3,4]) + (1-alpha)*sum(sqrt(sum(grad(2,:,:,:,:).^2,5)),[2,3,4]) + indicator(x);

end

