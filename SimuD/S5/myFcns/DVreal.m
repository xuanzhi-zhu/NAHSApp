function out = DVreal(xi,ref,parameters)
bv=parameters.bv;
bo=parameters.bo;
theta1=parameters.theta1;
theta2=parameters.theta2;
e1=parameters.e1;
e2=parameters.e2;
S=parameters.S;
M=parameters.M;
J=parameters.J;
D_l=parameters.D_l;
d_o=parameters.d_o;
d_q=parameters.d_q;

delta=parameters.delta;
k1=parameters.k1;
k2=parameters.k2;
k3=parameters.k3;
h1=parameters.h1;
h2=parameters.h2;
h3=parameters.h3;
c1=parameters.c1;
c2=parameters.c2;


pd0=ref(1:2);
pd1=ref(3:4);
pd2=ref(5:6);
pd3=ref(7:8);

p=xi(1:2);
v=xi(3:4);
q=xi(5:6);
o=xi(7);
hbv=xi(8:9);
hbo=xi(10);

invM=inv(M);
invJ=inv(J);
R=q2rot(q);

zp=p-pd0;
fp=R*v-pd1;
s=sigma(zp,parameters);
Ds=Dsigma(zp,parameters);
zv=v - R'*pd1 + k1.*invM*R'*s - delta;
phi=(M*S-S*M)*R'*pd1 - S*M*delta;
eta=invM*R'*s + k1.*R'*Ds*fp...
    + k2.*invM*zv - M*R'*pd2 - D_l*v...
    - d_q.*(e1'*v)*abs(e1'*v).*e1 + R'*hbv;
zo=e2'*(phi*o + eta);

tbv=hbv-bv;
tbo=hbo-bo;

out=[s', zv'*M*M, zo, h2.*tbv', h3.*tbo];

end

