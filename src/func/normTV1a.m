function norm = normTV1a(x,varargin)

if isempty(varargin)
    indicator = @(x) 0;
else
    indicator = varargin{1};
end

grad = L(x);
norm = sum(sqrt(sum(grad.^2,1)),[2,3,4,5]) + indicator(x);

end

