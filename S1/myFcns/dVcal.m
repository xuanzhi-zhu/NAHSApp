function out=dVcal(xi,ref_,parameters)
bv=parameters.bv;
bo=parameters.bo;
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
hp=parameters.hp;
hv=parameters.hv;
ho=parameters.ho;
% c1=parameters.c1;
% c2=parameters.c2;

pd0=ref_(1:2);
pd1=ref_(3:4);
pd2=ref_(5:6);
pd3=ref_(7:8);

p=xi(1:2);
v=xi(3:4);
q=xi(5:6);
o=xi(7);
hbv=xi(8:9);
hbo=xi(10);
huT=xi(11);
hutau=xi(12);

invM=inv(M);
invJ=inv(J);
R=q2rot(q);

zp=p-pd0;
fp=R*v-pd1;
s=sigma(zp,parameters);
Ds=Dsigma(zp,parameters);
DDs=DDsigma(zp,parameters);

zv=v - R'*pd1 + kp.*invM*R'*s - delta;

phi=(M*S-S*M)*R'*pd1 - S*M*delta;
eta=inv(hv).*invM*R'*s + kp.*R'*Ds*fp...
    + kv.*inv(hv).*invM*zv - M*R'*pd2 - D_l*v...
    - d_q.*(e1'*v)*abs(e1'*v).*e1 + R'*hbv;
zo=e2'*(phi*o + eta);

% mu=-kp.*s'*R*invM*R'*s...
%       - kv.*zv'*zv...
%       - ko.*zo*zo...
%       + s'*R*delta;

% out=mu;

dhbv=chiv(xi,ref_,parameters);
dhbo=chio(xi,ref_,parameters);

fp_prime=R*invM*((M*S-S*M)*v*o - D_l*v - d_q.*e1'*v*abs(e1'*v).*e1...
         + huT.*e1 - M*R'*pd2 + R'*hbv);

dphi=(M*S-S*M)*(-S*R'*pd1*o + R'*pd2);

fv_prime=-invM*S*M*(zv + R'*pd1 +delta)*o...
         + kp.*invM*R'*Ds*fp...
         + S*R'*pd1*o - R'*pd2 - invM*D_l*v...
         - d_q.*(e1'*v).*abs(e1'*v).*invM*e1...
         + huT.*invM*e1 + invM*R'*hbv;

eta_prime=R'*dhbv - S*R'*hbv*o...
          - inv(hv).*invM*S*R'*s*o + inv(hv).*invM*R'*Ds*fp...
          + kp.*R'*Ds*fp_prime - kp.*S*R'*Ds*fp*o...
          + kp.*R'*kron(fp',eye(2))*DDs*fp + kv.*inv(hv).*invM*fv_prime...
          + M*S*R'*pd2*o - M*R'*pd3 - (D_l + 2.*d_q.*abs(e1'*v).*e1*e1')...
          *invM*(-S*M*v*o + huT.*e1 - D_l*v - d_q.*e1'*v*abs(e1'*v).*e1 + R'*hbv);

fo_prime=e2'*(dphi*o+eta_prime) + invJ.*e2'*phi*(-d_o.*o+hbo+hutau);

mu=s'*fp...
     + hv.*zv'*M^2*fv_prime...
     + ho.*zo*fo_prime;

tbv=hbv-bv;
tbo=hbo-bo;
gammav=parameters.gammav;
gammao=parameters.gammao;

PetaPv=kp.*R'*Ds*R + kv.*inv(hv).*invM - D_l - 2.*d_q.*abs(e1'*v).*e1*e1';
muv=inv(gammav).*R* ( hv.*M*zv + ho.*invM*PetaPv'*e2*zo );

muo=inv(gammao).*ho.*invJ.*e2'*phi*zo;

out=mu+gammav*tbv'*(dhbv-muv)+gammao*tbo'*(dhbo-muo);

end
