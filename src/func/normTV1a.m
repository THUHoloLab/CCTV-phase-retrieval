function norm = normTV1a(x,varargin)

if isempty(varargin)
    indicator = @(x) 0;
else
    indicator = varargin{1};
end

grad = L(x);
norm = sum(sum(abs(grad),1),[2,3,4]) + indicator(x);

end

