clear all
close all

load data_
pd=ref(:,1:2);


layout = [1;1;2;2]*ones(1,8);
% h=create_axis(layout,18,...
%     'innerymargin',0.1,...
%     'botmargin',0,...
%     'innerxmargin',0.015,...
%     'leftmargin',0.1);
h=create_axis(layout,35,...
    'innerymargin',0.03,...
    'botmargin',0,...
    'innerxmargin',0.030,...%'innerxmargin',0.025,...
    'leftmargin',0.25);
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
aux=linspace(0,600,100);
plot(aux,bar_bv.*ones(length(aux),1),':','Color',lightgrey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(aux,bar_hbv.*ones(length(aux),1),':','Color',heavygrey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(aux,norm(bv).*ones(length(aux),1),'-','Color',red,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(t,(sum(hbv.^2,2)).^(0.5),'-','Color',black,'LineWidth',lineWidth/2,'MarkerSize',markerSize);hold off;
grid on

xlim([0,500])
ylim([0,6])

leg=legend({'$\bar{b}_v$','$\bar{\hat{b}}_v$','$|b_v|$','$|\hat{b}_v|$'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','northwest','NumColumns',1)
pos = leg.Position;
set(leg,...
    'Position',pos+[-0.28 -0.22 0 +0.25])

% xh=xlabel('','fontsize',12)
yh=ylabel('$[N]$','fontsize',12)

% Adjust position - negative values move left, positive move right
% Format: [x, y, z] in normalized units
% set(xh, 'Units', 'normalized');
% currentPos = get(xh, 'Position');
% set(xh, 'Position', [currentPos(1), currentPos(2)-0.05, currentPos(3)]);  % Move down

set(yh, 'Units', 'normalized');
currentPos = get(yh, 'Position');
set(yh, 'Position', [currentPos(1)-0.02, currentPos(2), currentPos(3)]);  % Move left
ax = gca; % current axes
ax.FontSize = 12;

%%%=================================
axes(h(2));
aux=linspace(0,600,100);
plot(aux,bar_bo.*ones(length(aux),1)*[-1 1],':','Color',lightgrey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(aux,bar_hbo.*ones(length(aux),1)*[-1 1],':','Color',heavygrey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(aux,bo.*ones(length(aux),1),'-','Color',yellow,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(t,abs(hbo),'-','Color',black,'LineWidth',lineWidth/2,'MarkerSize',markerSize);hold off;
grid on

xlim([0,500])
ylim([0,0.1])

leg=legend({'$\bar{b}_{\omega}$','','$\bar{\hat{b}}_{\omega}$','','$b_{\omega}$','$|\hat{b}_{\omega}|$'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','northwest')
pos = leg.Position;
set(leg,...
    'Position',pos+[-0.3 -0.22 +0 0.25])

xh=xlabel('$t\,[s]$','fontsize',12)
yh=ylabel('$[N]$','fontsize',12)

% Adjust position - negative values move left, positive move right
% Format: [x, y, z] in normalized units
set(xh, 'Units', 'normalized');
currentPos = get(xh, 'Position');
set(xh, 'Position', [currentPos(1), currentPos(2)-0.05, currentPos(3)]);  % Move down

set(yh, 'Units', 'normalized');
currentPos = get(yh, 'Position');
set(yh, 'Position', [currentPos(1)-0.05, currentPos(2), currentPos(3)]);  % Move left
ax = gca; % current axes
ax.FontSize = 12;


matlabfrag2('hatb_norm')