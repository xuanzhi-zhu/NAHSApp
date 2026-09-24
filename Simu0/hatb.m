clear all
close all

load data_
pd=ref(:,1:2);


layout = [1;1;2]*ones(1,4);
% h=create_axis(layout,18,...
%     'innerymargin',0.1,...
%     'botmargin',0,...
%     'innerxmargin',0.015,...
%     'leftmargin',0.1);
h=create_axis(layout,35,...
    'innerymargin',0.07,...
    'botmargin',0,...
    'innerxmargin',0.030,...%'innerxmargin',0.025,...
    'leftmargin',0.15);
colors = get(gca,'colororder');
% blue=colors(1,:);
blue=1/256.*[108 142 191];
% red=colors(2,:);
red=1/256.*[234 107 102];
lightred=1/256.*[194 107 102];
lightred=red./1.1;
% yellow=colors(3,:);
yellow=1/256.*[215 154 2];
purple=colors(4,:);
% green=colors(5,:);
green=1/256.*[130 179 102];
lightgrey=[0.7 0.7 0.7];
heavygrey=[0.3 0.3 0.3];
black=[0 0 0];
white=1/256.*[255 255 255];

lineWidth = 6;
markerSize = 20;

aux=linspace(0,2*pi,1000);
bar_circle_x=bar_hbv.*cos(aux);
bar_circle_y=bar_hbv.*sin(aux);

lbar_circle_x=bar_bv.*cos(aux);
lbar_circle_y=bar_bv.*sin(aux);

%%%=================================
axes(h(1));
plot(lbar_circle_x,lbar_circle_y,':','Color',lightgrey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(bar_circle_x,bar_circle_y,':','Color',heavygrey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(bv(1),bv(2),'square','MarkerEdgeColor',black,'MarkerFaceColor',red,'LineWidth',1,'MarkerSize',markerSize);hold on;
plot(hbv(:,1),hbv(:,2),'-','Color',black,'LineWidth',lineWidth/2,'MarkerSize',markerSize);hold on;
plot(hbv(1,1),hbv(1,2),'square','MarkerEdgeColor',black,'MarkerFaceColor',white,'LineWidth',1,'MarkerSize',markerSize);hold on;
plot(hbv(end,1),hbv(end,2),'hexagram','MarkerEdgeColor',black,'MarkerFaceColor',white,'LineWidth',1,'MarkerSize',markerSize);hold off;
axis equal
grid on


xlim([-6,6])
ylim([-6,6])

leg=legend({'$\bar{b}_v\ball$','$\bar{\hat{b}}_v\ball$','$b_v$','$\hat{b}_v$','$\hat{b}_v(0)$','$\hat{b}_v(600)$'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','best','NumColumns',1)
pos = leg.Position;
set(leg,...
    'Position',pos+[-0.65 -0.3 0 +0.35])

% xlabel('$x_I\,[m]$','fontsize',12)
% ylabel('$y_I\,[m]$','fontsize',12)
% ax = gca; % current axes
% ax.FontSize = 12;

%%%=================================
axes(h(2));
aux=linspace(0,600,100);
plot(aux,bar_bo.*ones(length(aux),1)*[-1 1],':','Color',lightgrey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(aux,bar_hbo.*ones(length(aux),1)*[-1 1],':','Color',heavygrey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(aux,bo.*ones(length(aux),1),'-','Color',yellow,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(t,hbo,'-','Color',black,'LineWidth',lineWidth/2,'MarkerSize',markerSize);hold off;
grid on

xlim([0,600])
ylim([-0.1,0.1])

leg=legend({'$\bar{b}_{\omega}\ball$','','$\bar{\hat{b}}_{\omega}\ball$','','$b_{\omega}$','$\hat{b}_{\omega}$'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','northeast')
pos = leg.Position;
set(leg,...
    'Position',pos+[0.03 +0.1 +0 0.25])

xlabel('$t\,[s]$','fontsize',12)
% ylabel('$$','fontsize',12)
ax = gca; % current axes
ax.FontSize = 12;


matlabfrag2('hatb')