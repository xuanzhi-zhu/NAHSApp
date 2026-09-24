function [dpsi,ref] = ref(psi,t,parameters)
c_max=parameters.c_max;
o1=parameters.o1;
psi1=parameters.psi1;
o2=parameters.o2;
psi2=parameters.psi2;
v_max=parameters.v_max;

% % % normal curve
% pd0=c_max.*[sin(o1*t+psi1);sin(o2*t+psi2)];
% pd1=c_max.*[o1*cos(o1*t+psi1);o2*cos(o2*t+psi2)];
% pd2=-diag([o1^2;o2^2])*pd0;
% pd3=-diag([o1^2;o2^2])*pd1;
% dpsi=0;%auxiliary

%Lissajous Curve
Dphi_=Dphi(psi,parameters);
DDphi_=DDphi(psi,parameters);
DDDphi_=DDDphi(psi,parameters);

dpsi=v_max*inv(norm(Dphi_));
ddpsi=-v_max.*(norm(Dphi_))^(-3).*Dphi_'*DDphi_*dpsi;

f=Drho(Dphi_);
g=DDphi_*dpsi;
h=DDphi_;
Df=DDrho(Dphi_)*DDphi_*dpsi;
Dh=DDDphi_*dpsi;
Dg=kron(dpsi,eye(2))*Dh + h*ddpsi;

pd0=c_max.*[sin(o1*psi+psi1);sin(o2*psi+psi2)];
pd1=v_max.*rho(Dphi_);
pd2=v_max.*Drho(Dphi_)*DDphi_*dpsi;
pd3=v_max.*( kron(g',eye(2))*Df + f*Dg );

ref=[pd0;pd1;pd2;pd3];

end

function out=Dphi(theta,parameters)
    c_max=parameters.c_max;
    o1=parameters.o1;
    psi1=parameters.psi1;
    o2=parameters.o2;
    psi2=parameters.psi2;
    out=c_max.*[o1*cos(o1*theta+psi1);o2*cos(o2*theta+psi2)];
end

function out=DDphi(theta,parameters)
    c_max=parameters.c_max;
    o1=parameters.o1;
    psi1=parameters.psi1;
    o2=parameters.o2;
    psi2=parameters.psi2;
    out=c_max.*[-o1*o1*sin(o1*theta+psi1);-o2*o2*sin(o2*theta+psi2)];
end

function out=DDDphi(theta,parameters)
    c_max=parameters.c_max;
    o1=parameters.o1;
    psi1=parameters.psi1;
    o2=parameters.o2;
    psi2=parameters.psi2;
    out=c_max.*[-o1*o1*o1*cos(o1*theta+psi1);-o2*o2*o2*cos(o2*theta+psi2)];
end

function out=rho(x)
    out=x./norm(x);
end

function out=Drho(x)
    out=inv(norm(x)).*(eye(2)-rho(x)*rho(x)');
end

function out=DDrho(x)
    rho_=rho(x);
    Drho_=Drho(x);
    out=-inv(norm(x)).*( reshape(Drho_,[4,1])*rho_' + ( kron(rho_,eye(2)) + kron(eye(2),rho_) )*Drho_ );
end