function n = normTV(x,lam)

g = D(x);
n = lam * norm1(g(:,:,1)) + lam * norm1(g(:,:,2));



function v = norm1(x)
    v = norm(x(:),1);
end

end
