function dpsi = F_ref(psi)
global parameters
v_max=parameters.v_max;

% % % normal curve
% pd0=c_max.*[sin(o1*t+psi1);sin(o2*t+psi2)];
% pd1=c_max.*[o1*cos(o1*t+psi1);o2*cos(o2*t+psi2)];
% pd2=-diag([o1^2;o2^2])*pd0;
% pd3=-diag([o1^2;o2^2])*pd1;
% dpsi=0;%auxiliary

%Lissajous Curve
Dphi_=Dphi(psi);

dpsi=v_max*inv(norm(Dphi_));

end

function out=Dphi(theta)
global parameters
    c_max=parameters.c_max;
    o1=parameters.o1;
    psi1=parameters.psi1;
    o2=parameters.o2;
    psi2=parameters.psi2;
    out=c_max.*[o1*cos(o1*theta+psi1);o2*cos(o2*theta+psi2)];
end