function utau=kappatau(xi,ref,parameters)

e1=parameters.e1;
e2=parameters.e2;
S=parameters.S;
M=parameters.M;
J=parameters.J;
D_l=parameters.D_l;
d_o=parameters.d_o;
d_q=parameters.d_q;

delta=parameters.delta;
kp=parameters.kp;
kv=parameters.kv;
ko=parameters.ko;
% hp=parameters.hp;
hv=parameters.hv;
ho=parameters.ho;


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
% invJ=inv(J);

uT=kappaT(xi,ref,parameters);

R=q2rot(q);
zp=p-pd0;
fp=R*v-pd1;
fp_prime=R*invM*((M*S-S*M)*v*o - D_l*v - d_q.*e1'*v*abs(e1'*v).*e1...
         + uT.*e1 - M*R'*pd2 + R'*hbv);
s=sigma(zp,parameters);
Ds=Dsigma(zp,parameters);
DDs=DDsigma(zp,parameters);

zv=v - R'*pd1 + kp.*invM*R'*s - delta;
fv_prime=-invM*S*M*(zv + R'*pd1 +delta)*o...
         + kp.*invM*R'*Ds*fp...
         + S*R'*pd1*o - R'*pd2 - invM*D_l*v...
         - d_q.*(e1'*v).*abs(e1'*v).*invM*e1...
         + uT.*invM*e1 + invM*R'*hbv;

phi=(M*S-S*M)*R'*pd1 - S*M*delta;
eta=inv(hv).*invM*R'*s + kp.*R'*Ds*fp...
    + kv.*inv(hv).*invM*zv - M*R'*pd2 - D_l*v...
    - d_q.*(e1'*v)*abs(e1'*v).*e1 + R'*hbv;
zo=e2'*(phi*o + eta);

dhbv=chiv(xi,ref,parameters);
dphi=(M*S-S*M)*(-S*R'*pd1*o + R'*pd2);
eta_prime=R'*dhbv - S*R'*hbv*o...
          - inv(hv).*invM*S*R'*s*o + inv(hv).*invM*R'*Ds*fp...
          + kp.*R'*Ds*fp_prime - kp.*S*R'*Ds*fp*o...
          + kp.*R'*kron(fp',eye(2))*DDs*fp + kv.*inv(hv).*invM*fv_prime...
          + M*S*R'*pd2*o - M*R'*pd3 - (D_l + 2.*d_q.*abs(e1'*v).*e1*e1')...
          *invM*(-S*M*v*o + uT.*e1 - D_l*v - d_q.*e1'*v*abs(e1'*v).*e1 + R'*hbv);

utau=d_o*o - hbo...
     - inv(ho).*J.*inv(e2'*phi).*e2'*(hv.*M*zv + ko.*zo.*e2 + ho.*dphi*o + ho.*eta_prime);
end

