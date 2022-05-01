function re = relative_error_2d(x_es,x_gt,region)

x_es = x_es(region.x1:region.x2,region.y1:region.y2);
x_gt = x_gt(region.x1:region.x2,region.y1:region.y2);

amp_es = abs(x_es);
pha_es = puma_ho(angle(x_es),1);

pha_gt = angle(x_gt);
pha_es = pha_es - mean(pha_es(:)) + mean(pha_gt(:));

x_es = amp_es.*exp(1i*pha_es);

re = norm2(x_es - x_gt) / norm2(x_gt);

function val = norm2(x)
    val = norm(x(:),2);
end

end

