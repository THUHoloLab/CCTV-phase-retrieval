function u = opr_transpose(x,params)

u = propagate(x,params.dist,params.pxsize,params.wavlen,params.method);
u = r2c_transpose(u);
u = real(u);

end

