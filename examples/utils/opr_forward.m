function u = opr_forward(x,params)

u = r2c(x);             % real vector -> complex vector
u = propagate(u,params.dist,params.pxsize,params.wavlen,params.method);
u = real(u);

end

