function [zp,zv,zo,tbv,tbo,Vreal] = errs(xi,ref,parameters)
e1=parameters.e1;
e2=parameters.e2;
S=parameters.S;
% m=parameters.m;
% mx=parameters.mx;
% my=parameters.my;
M=parameters.M;
D_l=parameters.D_l;
% d_o=parameters.d_o;
d_q=parameters.d_q;
% J=parameters.J;
bv=parameters.bv;
bo=parameters.bo;

delta=parameters.delta;
hp=parameters.hp;
hv=parameters.hv;
ho=parameters.ho;
kp=parameters.kp;
kv=parameters.kv;
ko=parameters.ko;
gammav=parameters.gammav;
gammao=parameters.gammao;

% c1=parameters.c1;

pd0=ref(1:2);
pd1=ref(3:4);
pd2=ref(5:6);
% pd3=ref(7:8);

p=xi(1:2);
v=xi(3:4);
q=xi(5:6);
o=xi(7);
hbv=xi(8:9);
hbo=xi(10);


invM=inv(M);
% invJ=inv(J);
R=q2rot(q);

zp=p-pd0;
fp=R*v-pd1;
s=sigma(zp,parameters);
Ds=Dsigma(zp,parameters);
zv=v - R'*pd1 + kp.*invM*R'*s - delta;
phi=(M*S-S*M)*R'*pd1 - S*M*delta;
eta=inv(hv).*invM*R'*s + kp.*R'*Ds*fp...
    + kv.*inv(hv).*invM*zv - M*R'*pd2 - D_l*v...
    - d_q.*(e1'*v)*abs(e1'*v).*e1 + R'*hbv;
zo=e2'*(phi*o + eta);

tbv=hbv-bv;
tbo=hbo-bo;



Vreal=hp.*(sqrt(norm(zp)^2+1) - 1)...
      + 0.5.*hv.*zv'*M*M*zv...
      + 0.5.*gammav.*tbv'*tbv...
      + 0.5.*ho.*zo*zo...
      + 0.5.*gammao.*tbo'*tbo;

% Vreal=max(0,Vreal);
end