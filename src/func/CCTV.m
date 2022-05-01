function val = CCTV(x,lam)

val = normTV(x,lam) + indicator(x);

end

