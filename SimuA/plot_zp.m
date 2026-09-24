clear all
close all

load data_S1 t zp b_prime
t_S1=t;
aux=size(zp);
zp_S1=reshape(zp,[aux(1) aux(3)]);
norm_zp_S1=(sum(zp_S1.^2,1)).^(0.5);
norm_zp_S1=norm_zp_S1';
clear t zp

load data_S2 t zp
t_S2=t;
aux=size(zp);
zp_S2=reshape(zp,[aux(1) aux(3)]);
norm_zp_S2=(sum(zp_S2.^2,1)).^(0.5);
norm_zp_S2=norm_zp_S2';
clear t zp

load data_S3 t zp
t_S3=t;
aux=size(zp);
zp_S3=reshape(zp,[aux(1) aux(3)]);
norm_zp_S3=(sum(zp_S3.^2,1)).^(0.5);
norm_zp_S3=norm_zp_S3';
clear t zp

load data_S4 t zp
t_S4=t;
aux=size(zp);
zp_S4=reshape(zp,[aux(1) aux(3)]);
norm_zp_S4=(sum(zp_S4.^2,1)).^(0.5);
norm_zp_S4=norm_zp_S4';
clear t zp

load data_S5 t zp
t_S5=t;
aux=size(zp);
zp_S5=reshape(zp,[aux(1) aux(3)]);
norm_zp_S5=(sum(zp_S5.^2,1)).^(0.5);
norm_zp_S5=norm_zp_S5';
clear t zp

load data_S6 t zp
t_S6=t;
aux=size(zp);
zp_S6=reshape(zp,[aux(1) aux(3)]);
norm_zp_S6=(sum(zp_S6.^2,1)).^(0.5);
norm_zp_S6=norm_zp_S6';
clear t zp

load data_S7 t zp
t_S7=t;
aux=size(zp);
zp_S7=reshape(zp,[aux(1) aux(3)]);
norm_zp_S7=(sum(zp_S7.^2,1)).^(0.5);
norm_zp_S7=norm_zp_S7';
clear t zp



%%%statistics
T_ss=30;
T_ss_fi=100;
Ind_S1=find(t_S1>T_ss, 1 );
Ind_S2=find(t_S2>T_ss, 1 );
Ind_S3=find(t_S3>T_ss, 1 );
Ind_S4=find(t_S4>T_ss, 1 );
Ind_S5=find(t_S5>T_ss, 1 );
Ind_S6=find(t_S6>T_ss, 1 );
Ind_S7=find(t_S7>T_ss, 1 );

Ind_S1_fi=find(t_S1<T_ss_fi, 1,'last' );
Ind_S2_fi=find(t_S2<T_ss_fi, 1,'last' );
Ind_S3_fi=find(t_S3<T_ss_fi, 1,'last' );
Ind_S4_fi=find(t_S4<T_ss_fi, 1,'last' );
Ind_S5_fi=find(t_S5<T_ss_fi, 1,'last' );
Ind_S6_fi=find(t_S6<T_ss_fi, 1,'last' );
Ind_S7_fi=find(t_S7<T_ss_fi, 1,'last' );

% t_S1=t_S1(Ind_S1:Ind_S1_fi);
% t_S2=t_S2(Ind_S2:Ind_S2_fi);
% t_S3=t_S3(Ind_S3:Ind_S1_fi);
% t_S4=t_S4(Ind_S4:Ind_S1_fi);
% t_S5=t_S5(Ind_S5:Ind_S1_fi);
% t_S6=t_S6(Ind_S6:Ind_S1_fi);
% t_S7=t_S7(Ind_S7:Ind_S1_fi);

norm_zp_S1_sta=norm_zp_S1(Ind_S1:Ind_S1_fi);
norm_zp_S2_sta=norm_zp_S2(Ind_S2:Ind_S2_fi);
norm_zp_S3_sta=norm_zp_S3(Ind_S3:Ind_S3_fi);
norm_zp_S4_sta=norm_zp_S4(Ind_S4:Ind_S4_fi);
norm_zp_S5_sta=norm_zp_S5(Ind_S5:Ind_S5_fi);
norm_zp_S6_sta=norm_zp_S6(Ind_S6:Ind_S6_fi);
norm_zp_S7_sta=norm_zp_S7(Ind_S7:Ind_S7_fi);


RMSE_S1=rms(norm_zp_S1_sta)
RMSE_S2=rms(norm_zp_S2_sta)
RMSE_S3=rms(norm_zp_S3_sta)
RMSE_S4=rms(norm_zp_S4_sta)
RMSE_S5=rms(norm_zp_S5_sta)
RMSE_S6=rms(norm_zp_S6_sta)
RMSE_S7=rms(norm_zp_S7_sta)

RMSE=[RMSE_S1;RMSE_S2;RMSE_S3;RMSE_S4;RMSE_S5;RMSE_S6;RMSE_S7]

%%%



% %find the index of an event
% indices=find((j-[0;j(1:end-1)])==1);
% %corresponding t value
% t_event=t(indices);
% %inter-event intervals
% t_inter=t_event-[0;t_event(1:end-1)];


% layout = [1;1;1;1]*ones(1,8);
% h=create_axis(layout,18,...
%     'innerymargin',0.1,...
%     'botmargin',0,...
%     'innerxmargin',0.08,...
%     'leftmargin',0.02);

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


lineWidth = 2;
markerSize = 3;


t_aux=linspace(0,600,10);
% zp_ub=b_prime^(0.5)*(1-b_prime)^(-0.5);
zp_ub=4.338;%a_brime

axes(h(1))

h7=plot(t_S7,norm_zp_S7,...
    'LineStyle',lines{1},'Color',yellow,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',yellow);hold on;
h6=plot(t_S6,norm_zp_S6,...
    'LineStyle',lines{1},'Color',red,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',red);hold on;
h5=plot(t_S5,norm_zp_S5,...
    'LineStyle',lines{1},'Color',purple,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',purple);hold on;
h4=plot(t_S4,norm_zp_S4,...
    'LineStyle',lines{1},'Color',green,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',green);hold on;
h3=plot(t_S3,norm_zp_S3,...
    'LineStyle',lines{1},'Color',blue,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',blue);hold on;
h2=plot(t_S2,norm_zp_S2,...
    'LineStyle',lines{1},'Color',grey,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',grey);hold on;
h1=plot(t_S1,norm_zp_S1,...
    'LineStyle',lines{1},'Color',black,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',black);hold on;
h0=plot(t_aux,zp_ub.*ones(numel(t_aux),1),...
    'LineStyle',lines{2},'Color',grey,'LineWidth',2*lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',grey);hold off;

grid on

leg=legend([h7 h6 h5 h4 h3 h2 h1 h0], {'S7)','S6)','S5)','S4)','S3)','S2)','S1)','$\sqrt{\bar{a}}$'},'fontsize',12);
% leg=legend([h1 h2 h3 h4 h5 h6 h7 h0], {'S1)','S2)','S3)','S4)','S5)','S6)','S7)','$\bar{a}$'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','best')
pos = leg.Position;
set(leg,...
    'Position',pos+[-0.83 -0.5 +0.05 0.35])

leg.Direction = 'reverse';


xlim([0,100])
ylim=[-5 70];
set(gca,'ylim',ylim,'YScale', 'linear')
% yticks([10^(-3) 10^(-2) 10^(-1) 1 10 100]);
% yticklabels({'0.001','0.01','0.1','1','10','100'})


% xh=xlabel('$t\,[s]$','fontsize',12)
yh=ylabel('$|p-p_d|\,[m]$','fontsize',12)

% Adjust position - negative values move left, positive move right
% Format: [x, y, z] in normalized units
% set(xh, 'Units', 'normalized');
% currentPos = get(xh, 'Position');
% set(xh, 'Position', [currentPos(1), currentPos(2)-0.05, currentPos(3)]);  % Move down

set(yh, 'Units', 'normalized');
currentPos = get(yh, 'Position');
set(yh, 'Position', [currentPos(1)-0.05, currentPos(2), currentPos(3)]);  % Move left
ax = gca; % current axes
ax.FontSize = 12;


axes(h(2))

lineWidth = 2;

h7=plot(t_S7,norm_zp_S7,...
    'LineStyle',lines{1},'Color',yellow,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',yellow);hold on;
h6=plot(t_S6,norm_zp_S6,...
    'LineStyle',lines{1},'Color',red,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',red);hold on;
h5=plot(t_S5,norm_zp_S5,...
    'LineStyle',lines{1},'Color',purple,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',purple);hold on;
h4=plot(t_S4,norm_zp_S4,...
    'LineStyle',lines{1},'Color',green,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',green);hold on;
h3=plot(t_S3,norm_zp_S3,...
    'LineStyle',lines{1},'Color',blue,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',blue);hold on;
h2=plot(t_S2,norm_zp_S2,...
    'LineStyle',lines{1},'Color',grey,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',grey);hold on;
h1=plot(t_S1,norm_zp_S1,...
    'LineStyle',lines{1},'Color',black,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',black);hold on;
h0=plot(t_aux,zp_ub.*ones(numel(t_aux),1),...
    'LineStyle',lines{2},'Color',grey,'LineWidth',2*lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',grey);hold off;

grid on
% 
% leg=legend([h7 h6 h5 h4 h3 h2 h1 h0], {'S7)','S6)','S5)','S4)','S3)','S2)','S1)','$\sqrt{\bar{a}}$'},'fontsize',12);
% % leg=legend([h1 h2 h3 h4 h5 h6 h7 h0], {'S1)','S2)','S3)','S4)','S5)','S6)','S7)','$\bar{a}$'},'fontsize',12);
% legend('boxoff')
% leg.Orientation='vertical';
% set(leg,...
%     'Location','best')
% pos = leg.Position;
% set(leg,...
%     'Position',pos+[-0.83 -0.5 +0.05 0.35])
% 
% leg.Direction = 'reverse';


xlim([0,100])
ylim=[-0.05 0.75];
set(gca,'ylim',ylim,'YScale', 'linear')
% yticks([10^(-3) 10^(-2) 10^(-1) 1 10 100]);
% yticklabels({'0.001','0.01','0.1','1','10','100'})


xh=xlabel('$t\,[s]$','fontsize',12)
yh=ylabel('$|p-p_d|\,[m]$','fontsize',12)

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

ax = gca; % current axes
ax.FontSize = 12;
matlabfrag2('simuA_zp')