function out = ln(x,bar_ln)

out=bar_ln.*log(1+inv(bar_ln).*x);

end

