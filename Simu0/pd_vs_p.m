clear all
close all

load data_
pd=ref(:,1:2);

ind_end=length(t);
% ind_end=find(t>=500,1);


layout = [1;1;1;1]*ones(1,8);
% h=create_axis(layout,15,...
%     'innerymargin',0.015,...
%     'botmargin',.13,...
%     'innerxmargin',0.06,...
%     'leftmargin',0.05);
h=create_axis(layout,18,...
    'innerymargin',0.1,...
    'botmargin',0,...
    'innerxmargin',0.015,...
    'leftmargin',0.1);
colors = get(gca,'colororder');
% blue=colors(1,:);
blue=1/256.*[108 142 191];
red=colors(2,:);
yellow=colors(3,:);
purple=colors(4,:);
% green=colors(5,:);
green=1/256.*[130 179 102];
grey=[0.7 0.7 0.7];
black=[0 0 0];

lineWidth = 1.5;
markerSize = 5;




plot(pd(1:ind_end,1),pd(1:ind_end,2),'--','Color',green,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(p(1:ind_end,1),p(1:ind_end,2),'-','Color',blue,'LineWidth',lineWidth/2,'MarkerSize',markerSize);hold on;
plot(pd(1,1),pd(1,2),'square','MarkerEdgeColor',black,'MarkerFaceColor',green,'LineWidth',0.5,'MarkerSize',markerSize);hold on;
plot(p(1,1),p(1,2),'square','MarkerEdgeColor',black,'MarkerFaceColor',blue,'LineWidth',0.5,'MarkerSize',markerSize);hold on;
plot(pd(ind_end,1),pd(ind_end,2),'hexagram','MarkerEdgeColor',black,'MarkerFaceColor',green,'LineWidth',0.5,'MarkerSize',markerSize);hold on;
plot(p(ind_end,1),p(ind_end,2),'hexagram','MarkerEdgeColor',black,'MarkerFaceColor',blue,'LineWidth',0.5,'MarkerSize',markerSize);hold off;
axis equal
grid on


xlim([-40,50])
ylim([-35,60])

leg=legend({'$p_d$','$p$','$p_d(0)$','$p(0)$','$p_d(500)$','$p(500)$'},'fontsize',12);
legend('boxoff')
leg.Orientation='horizontal';
set(leg,...
    'Location','best','NumColumns',2)
pos = leg.Position;
set(leg,...
    'Position',pos+[0.25 +0.08 -0.1 0])

xlabel('$x_I\,[m]$','fontsize',12)
ylabel('$y_I\,[m]$','fontsize',12)
ax = gca; % current axes
ax.FontSize = 12;

matlabfrag2('xy')