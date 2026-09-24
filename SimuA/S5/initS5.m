clear all
close all
clc

rng('shuffle');
addpath('./myFcns/')

ETC_on=0;

%Froude number<0.4 ===> |dpd|\leq 0.4*sqrt(9.81*L);
Fn=0.4;
L=1;%effective submerged length
% v_max=Fn*sqrt(9.81*L);

% %ASV (p,v,q,o)\in\reals[2]\x\reals[2]\x\sphere(1)\x\reals
% dp=R(q)*v;
% dv=-inv(M)*S*M*v*o - inv(M)*Dl*v - inv(M)*dq.*e1'*v*abs(e1'*v).*e1 + inv(M)*R(q)'*bv + inv(M)*uT*e1;
% dq=0.5.*S*q*o;
% do=-inv(J)*do*o + inv(J)*bo + inv(J)*utau;

%%%%perturbations (bv,bo)\in\reals[2]\x\reals
% bv=[10/3;20/3];
bv=0.5.*[10/3;20/3];
bo=0.05;
% bv=zeros(2,1);
% bo=0;
%%%%perturbations ub
bar_bv=1.25*norm(bv);
bar_bo=1.25*norm(bo);
%%%%estimation ub
bar_hbv=1.25*bar_bv;
bar_hbo=1.25*bar_bo;
%used in Projection operator
delta_bv=1e2;%the larger the smaller overshoot
delta_bo=1e3;

% %%%% ref
% % normal curve: pd0=c_max.*[sin(o1*t+psi1);sin(o2*t+psi2)];
% c_max=30;
% o1=1/60;
% psi1=pi/2;
% o2=1/60;
% psi2=0;
% v_max=1;%useless
% psi0=0;%useless
%Lissajous Curve: c_max.*[sin(o1*psi+psi1);sin(o2*psi+psi2)];
c_max=30;
o1=3;
psi1=pi/2;
o2=2;
psi2=0;
v_max=1;
psi0=0;

global parameters

% Parameters
parameters = struct();
parameters.ETC_on=ETC_on;
parameters.e1=[1;0];
parameters.e2=[0;1];
parameters.S=[0,-1;1,0];
parameters.m=25;
parameters.mx=-1.5;
parameters.my=-5.6;
parameters.M=diag([parameters.m+abs(parameters.mx);parameters.m+abs(parameters.my)]);
parameters.d_l=[25;37.5];
parameters.D_l=diag(parameters.d_l);
parameters.d_o=107;
parameters.d_q=7;
parameters.J=4.4;
parameters.bv=bv;
parameters.bo=bo;
parameters.bar_bv=bar_bv;
parameters.bar_bo=bar_bo;
parameters.bar_hbv=bar_hbv;
parameters.bar_hbo=bar_hbo;
parameters.delta_bv=delta_bv;
parameters.delta_bo=delta_bo;
parameters.eigs_M=eigs(parameters.M);

uThreshold=2;

parameters.uThreshold=uThreshold;

parameters.c_max=c_max;
parameters.o1=o1;
parameters.psi1=psi1;
parameters.o2=o2;
parameters.psi2=psi2;
parameters.v_max=v_max;

% parameters.delta=[-0.5;0];
% parameters.hp=1.5;%1.5 increase hp ===>smaller |zp|
% parameters.hv=0.5;%
% parameters.ho=0.5;%
% parameters.kp=30;%30
% parameters.kv=50;%20
% parameters.ko=10;%2
% 
% parameters.gammav=0.5;%increase gammav ===>slower estimation
% parameters.gammao=5;%increase gammao ===>slower estimation

parameters.delta=[-2;0];%[-0.5;0];
parameters.hp=4.5;%3
parameters.hv=0.01;%0.01
parameters.ho=0.15;%0.15
parameters.kp=30;%35
parameters.kv=10;%10
parameters.ko=5;%5

parameters.gammav=15;%5
parameters.gammao=1500;%1000

%check Ass on delta and Froude number, compactness
m=parameters.m;
mx=parameters.mx;
my=parameters.my;
delta=parameters.delta;
aux=min(abs(delta(1))*(m+abs(mx))/abs(mx-my),Fn*sqrt(9.81*L));
kp=parameters.kp;
kv=parameters.kv;
ko=parameters.ko;
hp=parameters.hp;
gammav=parameters.gammav;
gammao=parameters.gammao;
M=parameters.M;


% theta1=parameters.theta1;
% theta2=parameters.theta2;

% a_aux=inv(k1).*h1^(-2).*max(eigs(M));
% c_aux=a_aux.*( h1.*norm(delta) + theta1.*bar_bv.*bar_bv + theta2.*bar_bo.*bar_bo);
% c1_base=h1*(inv(sqrt(1-c_aux))-1)...
%      + ( 0.5.*max(eigs(M)).*max(eigs(M)).*inv(k2) + 0.5*inv(k3) + 2./3.*inv(theta1).*h2 + 2./3.*inv(theta2).*h3 ).*...
%      ( h1.*norm(delta) + theta1.*bar_bv.*bar_bv + theta2.*bar_bo.*bar_bo);
% 
% c1_scaling=3;
% c1=c1_scaling*c1_base;
% parameters.c1=c1;
% 
% c2=0.5;%\in (0,1)
% parameters.c2=c2;


% %normal curve
% flag1=c_max*(o1^2+o2^2)^(0.5)-aux;
%Lissajous Curve
flag1=v_max-aux;

flag2=kp-hp^(-2)*max(eigs(M))*( hp*norm(delta) );

if flag1>=0
    % disp('decrease ref velocities');
    error('decrease ref velocities');
end

if flag2<=0
    % disp('decrease ref velocities');
    error('increase k1');
end

% Initial conditions
p0=[45;20];
v0=[0;0];
ang0=pi;%heading angle of {B} wrt {I}
R0=ang2rot(ang0);
q0=ang2q(ang0);
q0=q0./norm(q0);
o0=0;
hbv0=[0;0];%inside the ball of radius bar_hbv
hbo0=0;%inside the ball of radius bar_hbo

if (norm(hbv0)>=bar_hbv)||(norm(hbo0)>=bar_hbo)
    % disp('decrease ref velocities');
    error('decrease hbv0 or hbo0');
end

%plant state
x0 = [p0;v0;q0;o0;hbv0;hbo0];

timer0=0;

%ref
[~,ref0]=ref(psi0,timer0,parameters);
%memory state (unused here)
huT0=0;
hutau0=0;

%check Lya (V=max{0,V_z-c1})
[~,~,~,~,~,Vcal0]=errs(x0,ref0,parameters);


xi0=[x0;huT0;hutau0;timer0;Vcal0];


% Simulation horizon
TSPAN = [0 620];
JSPAN = [0 30*TSPAN(2)];
% rule=2;%2 for flow
                                                                        
% Solver tolerances
AbsTol = 1e-3;
RelTol = AbsTol*AbsTol;
MaxStep = 1e-2;%

% options=odeset('AbsTol',AbsTol,'RelTol',RelTol,'MaxStep',MaxStep);
% 
% [t,j,xi]=HyEQsolver(@F,@G,@C,@D,xi0,TSPAN,JSPAN,rule,options);
tic 
sim('modelS5.slx')
toc

save data_S5