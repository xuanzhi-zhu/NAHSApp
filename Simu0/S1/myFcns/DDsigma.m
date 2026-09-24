function out=DDsigma(zp,parameters)
    hp=parameters.hp;
    s=sigma(zp,parameters);
    Ds=Dsigma(zp,parameters);
    out= - hp^(-1)./sqrt(norm(zp)^2+1)...
        .* (reshape(Ds,[4,1])*s' + (kron(s,eye(2))+kron(eye(2),s))*Ds);
end
