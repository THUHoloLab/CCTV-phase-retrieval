function y = proj(x)

x = r2c(x);
arg = angle(x);
mod = abs(x);
y = 1.0*mod./max(1.0,mod).*exp(1i*arg);
y = c2r(y);

end

