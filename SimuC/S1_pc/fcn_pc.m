function pc=fcn_pc(scale_i,pc_nom,parameters)

m=parameters.m;
mx=parameters.mx;
my=parameters.my;
delta=parameters.delta;
d_l=parameters.d_l;
M=parameters.M;
J=parameters.J;
e1=parameters.e1;
e2=parameters.e2;
bar_bv=parameters.bar_bv;
bar_bo=parameters.bar_bo;
bar_hbv=parameters.bar_hbv;
bar_hbo=parameters.bar_hbo;
v_max=parameters.v_max;

s=scale_i+1;

hp_nom=pc_nom(1);
hv_nom=pc_nom(2);
ho_nom=pc_nom(3);
kp_nom=pc_nom(4);
kv_nom=pc_nom(5);
ko_nom=pc_nom(6);

hp=hp_nom*s;
hv=hv_nom;
ho=1/(kp_nom*hp_nom*s);
kp=kp_nom;
kv=kv_nom*s;

delta1=ho*(min(eigs(M)))^(-1)*(kp*hp+e2'*d_l)*(bar_bv+bar_hbv) ...
             + ho*J^(-1)*((m+abs(mx))*abs(delta(1))+abs(mx-my)*v_max)*(bar_bo+bar_hbo);
delta2=ho*(min(eigs(M)))^(-1)*hv*(min(eigs(M)))^(-1)*(bar_bv+bar_hbv);

ko=s*(delta1+kv_nom*s*delta2)^2;

pc=[hp;hv;ho;kp;kv;ko];

end

