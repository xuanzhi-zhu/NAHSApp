function out=sigma(zp,parameters)
    hp=parameters.hp;
    out= hp./sqrt(norm(zp)^2+1).*zp;
end
