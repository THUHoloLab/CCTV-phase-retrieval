function val = normTV2i(x,alpha,varargin)

if isempty(varargin)
    indicator = @(x) 0;
else
    indicator = varargin{1};
end

grad = L(x);
val = alpha*sum(sqrt(sum(grad(1:3,:,:,:).^2,1)) + (1-alpha)*sqrt(sum(grad(4:6,:,:,:).^2,1)),[2,3,4]) + indicator(x);

end

