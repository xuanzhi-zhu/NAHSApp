function [c,ceq] = cons(z)
global opti_para
delta=opti_para.delta;
k1=opti_para.k1;
k2=opti_para.k2;
k3=opti_para.k3;
h1=opti_para.h1;
M=opti_para.M;
theta1=opti_para.theta1;
theta2=opti_para.theta2;
bar_bv=opti_para.bar_bv;
bar_bo=opti_para.bar_bo;

c = k1.*h1.*h1.*inv(max(eigs(M))).*norm(z(1:2))^2./(norm(z(1:2))^2+1)...
                    + k2.*norm(z(3:4))^2 + k3.*z(7)^2 + (theta1-0.5.*z(9)).*norm(z(5:6))^2 + (theta2-0.5.*z(10)).*z(8)*z(8)...
                    - h1.*norm(delta) - 0.5.*inv(z(9)).*theta1^2.*bar_bv^2 - 0.5.*inv(z(10)).*theta2^2.*bar_bo^2;
ceq = [];
end

