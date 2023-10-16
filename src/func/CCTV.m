function val = CCTV(x,lam)
% =========================================================================
% The constrained complex total variation (CCTV) function.
% -------------------------------------------------------------------------
% Input:    - x   : The complex-valued 2D transmittance of the sample.
%           - lam : The regularization parameter for the total variation
%                   function.
% Output:   - val : Value of the CCTV function.
% =========================================================================

val = normTV(x,lam) + indicator(x);

end

