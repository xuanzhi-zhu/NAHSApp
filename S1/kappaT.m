function uT=kappaT(xi,ref,parameters)
e1=parameters.e1;
% e2=parameters.e2;
S=parameters.S;
M=parameters.M;
J=parameters.J;
D_l=parameters.D_l;
% d_o=parameters.d_o;
d_q=parameters.d_q;

delta=parameters.delta;
kp=parameters.kp;
kv=parameters.kv;
% ko=parameters.ko;
% hp=parameters.hp;
hv=parameters.hv;
% ho=parameters.ho;

pd0=ref(1:2);
pd1=ref(3:4);
pd2=ref(5:6);
% pd3=ref(7:8);

p=xi(1:2);
v=xi(3:4);
q=xi(5:6);
o=xi(7);
hbv=xi(8:9);
% hbo=xi(10);

invM=inv(M);
% invJ=inv(J);


R=q2rot(q);
zp=p-pd0;
fp=R*v-pd1;
s=sigma(zp,parameters);
Ds=Dsigma(zp,parameters);
zv=v - R'*pd1 + kp.*inv(M)*R'*s - delta;

phi=(M*S-S*M)*R'*pd1 - S*M*delta;

eta=inv(hv).*invM*R'*s + kp.*R'*Ds*fp...
    + kv.*inv(hv).*invM*zv - M*R'*pd2 - D_l*v...
    - d_q.*(e1'*v)*abs(e1'*v).*e1 + R'*hbv;

uT=-e1'*(phi*o + eta);
end