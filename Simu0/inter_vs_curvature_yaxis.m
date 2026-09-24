clear all
close all

load data_

%find the index of an event
indices=find((j-[0;j(1:end-1)])==1);
%corresponding t value
t_event=t(indices);
%inter-event intervals
t_inter=t_event-[0;t_event(1:end-1)];

% 
% PI=zeros(size(t_event));
% for i=1:1:numel(PI)
%     PI(i)=1/i*sum(t_inter(1:i));
% end


S=[0,-1;1,0];
%curvature
pd1=ref(:,3:4);
pd2=ref(:,5:6);
curvature=zeros(size(t));
for i=1:1:length(t)
    curvature(i)=(norm(pd1(i,:)))^(-3) * norm(pd2(i,:)*S*pd1(i,:)');
end


layout = [1;1;1;1]*ones(1,8);
h=create_axis(layout,18,...
    'innerymargin',0.1,...
    'botmargin',0,...
    'innerxmargin',0.08,...
    'leftmargin',0.02);
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

lineWidth = 3;
markerSize = 5;

plot(t,25.*curvature,'-','Color',grey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
plot(t_event(2:end),t_inter(2:end),'diamond','MarkerEdgeColor',black,'LineWidth',0.1,'MarkerSize',markerSize);hold off;

grid on
xh=xlabel('$t\,[s]$','fontsize',12)
yh=ylabel('$[s]$','fontsize',12)

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

xlim([0,500])
ylim([0.001,10])
set(gca,'ylim',ylim,'YScale', 'linear')

leg=legend({'$25\kappa$','$t_{j+1}-t_j$'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','northeast')
pos = leg.Position;
set(leg,...
    'Position',pos+[0.0 -0.05 +0 0.05])


ax = gca; % current axes
ax.FontSize = 12;
matlabfrag2('inter_curvature')