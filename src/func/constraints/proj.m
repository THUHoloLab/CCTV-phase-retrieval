function y = proj(x)

global constraint
global absorption
global support

if strcmpi(constraint,'none')
    y = x;
elseif strcmpi(constraint,'a')
    y = min(abs(x),absorption).*exp(1i*angle(x));
elseif strcmpi(constraint,'s')
    y = x.*support;
elseif strcmpi(constraint,'as')
    y = min(abs(x),absorption).*exp(1i*angle(x)).*support;
else
    error("Invalid constraint. Should be 'A'(absorption), 'S'(support), 'AS'(both), or 'none'.")
end

