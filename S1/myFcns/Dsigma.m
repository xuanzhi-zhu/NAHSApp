function out=Dsigma(zp,parameters)
    hp=parameters.hp;
    s=sigma(zp,parameters);
    out= hp./sqrt(norm(zp)^2+1).* (eye(2)-hp^(-2).*s*s');
end
