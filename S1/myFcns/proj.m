function out=proj(mu,hb, bar_b,bar_hb,delta)

out=mu-eta1(hb, bar_b,bar_hb)*eta2(mu,hb, bar_b,delta)*inv(bar_b).*hb;

end

function out=eta1(s, bar_b,bar_hb)

dist=max(0,norm(s)-bar_b);
out=(bar_hb-bar_b)^(-2)*dist^2;

end


function out=eta2(t,s, bar_b,delta)

aux=inv(bar_b)*t'*s;
out=aux+ (aux^2+4*delta^2)^(0.5);
out=0.5*out;

end
