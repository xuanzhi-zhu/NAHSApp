load data_S1
addpath('./myFcns/')

% p=xi(:,1:2);
% v=xi(:,3:4);
% q=xi(:,5:6);
% o=xi(:,7);
% hbv=xi(:,8:9);
% hbo=xi(:,10);
% huT=xi(:,11);
% hutau=xi(:,12);
% timer=xi(:,13);
% psi=xi(:,14);
% Vcal=xi(:,15);
% 
% pd0=zeros(length(t),2);
% 
% zp=zeros(length(t),2);
% zv=zeros(length(t),2);
% zo=zeros(length(t),1);
% tbv=zeros(length(t),2);
% tbo=zeros(length(t),1);
% Vreal=zeros(length(t),1);
% uT=zeros(length(t),1);
% utau=zeros(length(t),1);
% 
% for i=1:1:length(t)
%     [~,ref_] = ref(psi(i),timer(i),parameters);
%     pd0(i,:)=ref_(1:2)';
% 
%     [zp_,zv_,zo_,tbv_,tbo_,Vreal_] = errs(xi(i,:)',ref_,parameters);
%     zp(i,:)=zp_';
%     zv(i,:)=zv_';
%     zo(i,:)=zo_;
%     tbv(i,:)=tbv_';
%     tbo(i,:)=tbo_;
%     Vreal(i,:)=Vreal_;
%     uT(i)=kappaT(xi(i,:)',ref_,parameters);
%     utau(i)=kappatau(xi(i,:)',ref_,parameters);
% end

% %%%%t,j
% figure(1);plot(t,j)



% p=reshape(p,[2 length(tout)])';
% p=reshape(p,[2 length(tout)])';
% p=reshape(p,[2 length(tout)])';

% %%%%%errs
% figure(2);
% plot(t,zp);grid on
% legend('zp')
% 
% figure(3);
% plot(t,zv);grid on
% legend('zv')
% 
% figure(4);
% plot(t,zo);grid on
% legend('zo')
% 
% figure(5);
% plot(t,tbv);grid on
% legend('tbv')
% 
% figure(6);
% plot(t,tbo);grid on
% legend('tbo')

% %%%%Lya
% figure(7);
% plot(t,Vreal);hold on
% plot(t,Vcal);hold off
% legend('Vreal','Vcal')
% grid on;

% %%%%Lya
% figure(7);
% plot(t,abs(Vreal-Vcal));
% legend('Vreal-Vcal')
% set(gca, 'YScale', 'log')
% grid on;


% %%%controls
% figure(8);
% plot(t,huT);hold on
% plot(t,hutau);hold off
% legend('huT','hutau')
% grid on;


pd0=ref(:,1:2);
%%%%animation
figure(100);
Speedy=0.1/MaxStep;
for i=1:Speedy:length(t)
    %past locations
    plot(pd0(1:i,1),pd0(1:i,2),'ro','MarkerSize',2);hold on;
    plot(p(1:i,1),p(1:i,2),'bo','MarkerSize',1);hold on;
    %current location
    plot(pd0(i,1),pd0(i,2),'ro','MarkerSize',14);hold on;
    plot(p(i,1),p(i,2),'bo','MarkerSize',7);hold on;
    %initial location
    plot(pd0(1,1),pd0(1,2),'rx','MarkerSize',14);hold on;
    plot(p(1,1),p(1,2),'bx','MarkerSize',7);hold on;
    %current xB axis (e1) to be expressed in I
    R=q2rot(q(i,:));
    xB_I=10.*R*[1;0];
    arrow_ini=p(i,1:2)';
    arrow_end=arrow_ini + xB_I;
    plot(linspace(arrow_ini(1),arrow_end(1),10),linspace(arrow_ini(2),arrow_end(2),10),'LineWidth',2);hold on;
    %force disturbance bv
    arrow_ini=[0;0];
    arrow_end=arrow_ini + 10.*bv./norm(bv);
    plot(linspace(arrow_ini(1),arrow_end(1),10),linspace(arrow_ini(2),arrow_end(2),10),'LineWidth',2);hold on;
    %torque disturbance bo
    plot(p(i,1),p(i,2),'bx','MarkerSize',3); hold off;
    insidePaper=-abs(bo);
    text(p(i,1)+5,p(i,2)+5,num2str(insidePaper));

    grid on;axis equal;
    xlabel('x');ylabel('y');
    axis([-50 50 -50 50]);
    timer=tout(i);
    text(49,49,num2str(timer));
    %set(gca,'Ydir','reverse')
    pause(0.001)
    % pause(0.1)
end