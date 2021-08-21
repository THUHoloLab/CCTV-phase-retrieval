function val = indicator(x)

x = r2c(x);

if sum(sum(abs(x)>1.0+eps)) == 0
    val = 0;
else
    val = inf;
end

end

