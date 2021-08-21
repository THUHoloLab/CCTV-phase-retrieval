function grad = dF(x,y,params)

u = r2c(x);
u = propagate(u,params.dist,params.pxsize,params.wavlen,params.method);
r = abs(u) - sqrt(y);
u = r.*exp(1i*angle(u));
u = propagate(u,-params.dist,params.pxsize,params.wavlen,params.method);
u = r2c_hermitian(u);
grad = real(u);

end

