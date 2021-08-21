function norm = normTV2a(x,alpha,varargin)

if isempty(varargin)
    indicator = @(x) 0;
else
    indicator = varargin{1};
end

grad = L(x);
norm = alpha*sum(sum(abs(grad(1:3,:,:,:)),1) + (1-alpha)*sum(abs(grad(4:6,:,:,:)),1),[2,3,4]) + indicator(x);

end

