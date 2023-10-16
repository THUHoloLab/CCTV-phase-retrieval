function n = normTV(x,lam)
% =========================================================================
% Calculate the total variation (TV) seminorm for an input image.
% -------------------------------------------------------------------------
% Input:    - x   : 2D image.
%           - lam : Regularization weight.
% Output:   - n   : The funciton value.
% =========================================================================

g = D(x);
n = lam * norm1(g(:,:,1)) + lam * norm1(g(:,:,2));

function v = norm1(x)
    v = norm(x(:),1);
end

end
