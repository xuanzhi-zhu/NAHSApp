function [inC,inD]=fcn_decide(xi, ref, parameters)

e1=parameters.e1;
e2=parameters.e2;
S=parameters.S;
M=parameters.M;
J=parameters.J;
D_l=parameters.D_l;
d_o=parameters.d_o;
d_q=parameters.d_q;
eigs_M=parameters.eigs_M;

delta=parameters.delta;
kp=parameters.kp;
kv=parameters.kv;
ko=parameters.ko;
hp=parameters.hp;
hv=parameters.hv;
ho=parameters.ho;
% gammav=parameters.gammav;
% gammao=parameters.gammao;
c1=parameters.c1;
c2=parameters.c2;
bar_bv=parameters.bar_bv;
bar_hbv=parameters.bar_hbv;
bar_bo=parameters.bar_bo;
bar_hbo=parameters.bar_hbo;

d_l=parameters.d_l;
m=parameters.m;
mx=parameters.mx;
my=parameters.my;
v_max=parameters.v_max;

pd3_ub=parameters.pd3_ub;

% timer=xi(13);
% psi=xi(14);
% [~,ref_] = ref(psi,timer,parameters);

pd0=ref(1:2);
pd1=ref(3:4);
pd2=ref(5:6);
% pd3=ref_(7:8);

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

%===========================================================
Vz=hp.*(sqrt(norm(zp)^2+1) - 1)...
      + 0.5.*hv.*zv'*M*M*zv...
      + 0.5.*ho.*zo*zo;
%===========================================================

dhbv=chiv(xi,ref,parameters);
% dhbo=chio(xi,ref_,parameters);
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
          + M*S*R'*pd2*o - (D_l + 2.*d_q.*abs(e1'*v).*e1*e1')...
          *invM*(-S*M*v*o + huT.*e1 - D_l*v - d_q.*e1'*v*abs(e1'*v).*e1 + R'*hbv);
fo_prime=e2'*(dphi*o+eta_prime) + invJ.*e2'*phi*(-d_o.*o+hbo+hutau);

%===========================================================
DVz=[s', hv.*zv'*M^2, ho.*zo];
Fz=[fp;fv_prime;fo_prime];
%===========================================================

aux=norm(zv)*hv*norm(M)*(bar_bv+bar_hbv) ...
     + norm(zo)*ho*inv(min(eigs_M))*(kp*hp + kv*inv(hv)*inv(min(eigs_M)) + d_l(2))*(bar_bv+bar_hbv) ...
     + norm(zo)*ho*invJ*( (m+abs(mx))*abs(delta(1)) + abs(mx-my)*v_max)*(bar_bo+bar_hbo) ...
     - norm(ho*zo*e2'*M)*pd3_ub;

mu=-kp.*s'*R*invM*R'*s...
      - kv.*zv'*zv...
      - ko.*zo*zo...
      + s'*R*delta;

if (Vz<=c1) || (DVz*Fz<=c2.*mu + aux)
    inC=1;
else
    inC=0;
end

if (Vz>=c1) && (DVz*Fz >=c2.*mu + aux)
    inD=1;
else
    inD=0;
end

end