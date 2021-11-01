function val = normTV2i(x,alpha,varargin)

if isempty(varargin)
    indicator = @(x) 0;
else
    indicator = varargin{1};
end

grad = L(x);
val = alpha*sum(sqrt(sum(grad(1,:,:,:,:).^2,[2,3,4,5]))) + (1-alpha)*sum(sqrt(sum(grad(2,:,:,:,:).^2,[2,3,4,5]))) + indicator(x); 

end

