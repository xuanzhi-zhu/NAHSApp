clear all
close all

load data_S1 t j
t_S1=t;
j_S1=j;
%find the index of an event
indices=find((j_S1-[0;j_S1(1:end-1)])==1);
%corresponding t value
t_event_S1=t_S1(indices);
%inter-event intervals
t_inter_S1=t_event_S1-[0;t_event_S1(1:end-1)];
clear t j

load data_S2 t j
t_S2=t;
j_S2=j;
%find the index of an event
indices=find((j_S2-[0;j_S2(1:end-1)])==1);
%corresponding t value
t_event_S2=t_S2(indices);
%inter-event intervals
t_inter_S2=t_event_S2-[0;t_event_S2(1:end-1)];
clear t j

load data_S3 t j
t_S3=t;
j_S3=j;
%find the index of an event
indices=find((j_S3-[0;j_S3(1:end-1)])==1);
%corresponding t value
t_event_S3=t_S3(indices);
%inter-event intervals
t_inter_S3=t_event_S3-[0;t_event_S3(1:end-1)];
clear t j

load data_S4 t j
t_S4=t;
j_S4=j;
%find the index of an event
indices=find((j_S4-[0;j_S4(1:end-1)])==1);
%corresponding t value
t_event_S4=t_S4(indices);
%inter-event intervals
t_inter_S4=t_event_S4-[0;t_event_S4(1:end-1)];
clear t j

load data_S5 t j
t_S5=t;
j_S5=j;
%find the index of an event
indices=find((j_S5-[0;j_S5(1:end-1)])==1);
%corresponding t value
t_event_S5=t_S5(indices);
%inter-event intervals
t_inter_S5=t_event_S5-[0;t_event_S5(1:end-1)];
clear t j

load data_S6 t j
t_S6=t;
j_S6=j;
%find the index of an event
indices=find((j_S6-[0;j_S6(1:end-1)])==1);
%corresponding t value
t_event_S6=t_S6(indices);
%inter-event intervals
t_inter_S6=t_event_S6-[0;t_event_S6(1:end-1)];
clear t j

load data_S7 t j
t_S7=t;
j_S7=j;
%find the index of an event
indices=find((j_S7-[0;j_S7(1:end-1)])==1);
%corresponding t value
t_event_S7=t_S7(indices);
%inter-event intervals
t_inter_S7=t_event_S7-[0;t_event_S7(1:end-1)];
clear t j

layout = [[1;1]*ones(1,3) [2;2]*ones(1,3)];
h=create_axis(layout,13.4,...
    'innerymargin',0.05,...
    'botmargin',0,...
    'innerxmargin',0.1,...
    'leftmargin',0.03);
% colors = get(gca,'colororder');
blue=1/256.*[108 142 191];
green=1/256.*[130 179 102];

purple=1/256.*[189,119,149];

red=1/256.*[236,110,102];
yellow=1/256.*[247,172,83];

grey=[0.7 0.7 0.7];
black=[0 0 0];



lines={'-','--','-.',':'};% lines{1};
markers={'+','o','*','.','x','s','d','^','v','>','<','p','h','|','_'};% markers{1};


lineWidth = 0.6;
markerSize = 5;


%===========================================================

% aux_0=10;
% X1=linspace(-1.5,1,10);
% Y1=linspace(aux_0-6,aux_0-2,10);
% Y2=linspace(aux_0-8,aux_0-4,10);
% t_event_S2=linspace(0,30,10);
% t_inter_S2=aux_0.*ones(numel(t_event_S2),1);

axes(h(1))
h_aux7=plot(t_event_S7(2:end),t_inter_S7(2:end),...
    markers{13},'MarkerEdgeColor',yellow,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_aux6=plot(t_event_S6(2:end),t_inter_S6(2:end),...
    markers{9},'MarkerEdgeColor',red,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_aux5=plot(t_event_S5(2:end),t_inter_S5(2:end),...
    markers{8},'MarkerEdgeColor',purple,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_aux4=plot(t_event_S4(2:end),t_inter_S4(2:end),...
    markers{7},'MarkerEdgeColor',green,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_aux3=plot(t_event_S3(2:end),t_inter_S3(2:end),...
    markers{2},'MarkerEdgeColor',blue,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
% h_aux2=plot(t_event_S2,t_inter_S2,...
%     lines{1},'LineWidth',5*lineWidth,'Color',grey);hold on;
h_aux1=plot(t_event_S1(2:end),t_inter_S1(2:end),...
    markers{6},'MarkerEdgeColor',black,'LineWidth',lineWidth,'MarkerSize',markerSize);hold off;
% plot(X1,Y1,'-','Color',black,'lineWidth',3*lineWidth);hold on;
% plot(X1,Y2,'-','Color',black,'lineWidth',3*lineWidth);hold off;
% text(-5.5,8.5,'$+\infty$','fontsize',12)
grid on

xlim([0,20])
ylim([0.0001,10])

set(gca,'YScale', 'log')

leg=legend({'S7)','S6)','S5)','S4)','S3)','S1)'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','best','NumColumns',1)
pos = leg.Position;
set(leg,...
    'Position',pos+[0.15 0 0 0.05])


% ytickLabels = cellstr(num2str(round(log10(yticks(:))), '10^%d'));

yticks = ([10^(-4) 10^(-2) 10^(0)]);
yticklabels({'10^{-4}','10^{-2}','10^{0}'})



leg.Direction = 'reverse';
leg.ItemTokenSize = [5,5]; 

xlabel('$t\,[s]$','fontsize',12);
ylabel('$t_{j+1}-t_j\,[s]$','fontsize',12)





% aux_0=10;
% X1=linspace(580-1.5,580+1,10);
% Y1=linspace(aux_0-4,aux_0-1.5,10);
% Y2=linspace(aux_0-6,aux_0-3,10);
% t_event_S2=linspace(580,600,10);
% t_inter_S2=aux_0.*ones(numel(t_event_S2),1);

axes(h(2))
h_aux7=plot(t_event_S7(2:end),t_inter_S7(2:end),...
    markers{13},'MarkerEdgeColor',yellow,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_aux6=plot(t_event_S6(2:end),t_inter_S6(2:end),...
    markers{9},'MarkerEdgeColor',red,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_aux5=plot(t_event_S5(2:end),t_inter_S5(2:end),...
    markers{8},'MarkerEdgeColor',purple,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_aux4=plot(t_event_S4(2:end),t_inter_S4(2:end),...
    markers{7},'MarkerEdgeColor',green,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_aux3=plot(t_event_S3(2:end),t_inter_S3(2:end),...
    markers{2},'MarkerEdgeColor',blue,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
% h_aux2=plot(t_event_S2,t_inter_S2,...
%     lines{1},'LineWidth',5*lineWidth,'Color',grey);hold on;
h_aux1=plot(t_event_S1(2:end),t_inter_S1(2:end),...
    markers{6},'MarkerEdgeColor',black,'LineWidth',lineWidth,'MarkerSize',markerSize);hold off;
% plot(X1,Y1,'-','Color',black,'lineWidth',3*lineWidth);hold on;
% plot(X1,Y2,'-','Color',black,'lineWidth',3*lineWidth);hold off;
% text(574.5,8.5,'$+\infty$','fontsize',12)

grid on
xlim([480,500])
ylim([0.01,10])
set(gca,'ylim',ylim,'YScale', 'log')


yticks = ([10^(-2) 10^(-1) 10^(0)]);
yticklabels({'10^{-2}','10^{-1}','10^{0}'})

xlabel('$t\,[s]$','fontsize',12);



ax = gca; % current axes
ax.FontSize = 10;
matlabfrag2('simuA_inter')