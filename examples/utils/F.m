function val = F(x,y,params)

u = r2c(x);
u = propagate(u,params.dist,params.pxsize,params.wavlen,params.method);
val = 0.5*normArr(abs(u)-sqrt(y))^2;

end

