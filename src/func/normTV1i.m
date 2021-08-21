function val = normTV1i(x,varargin)

if isempty(varargin)
    indicator = @(x) 0;
else
    indicator = varargin{1};
end

grad = L(x);
val = sum(sqrt(sum(grad.^2,1)),[2,3,4]) + indicator(x);

end

