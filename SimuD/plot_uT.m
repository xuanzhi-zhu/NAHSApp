clear all
close all

load data_S1 t huT hutau
t_S1=t;
huT_S1=huT;
hutau_S1=hutau;
clear t huT hutau

load data_S2 t uT utau
t_S2=t;
huT_S2=uT;
hutau_S2=utau;
clear t uT utau

load data_S3 t huT hutau
t_S3=t;
huT_S3=huT;
hutau_S3=hutau;
clear t huT hutau

load data_S4 t huT hutau
t_S4=t;
huT_S4=huT;
hutau_S4=hutau;
clear t huT hutau

load data_S5 t huT hutau
t_S5=t;
huT_S5=huT;
hutau_S5=hutau;
clear t huT hutau

load data_S6 t huT hutau
t_S6=t;
huT_S6=huT;
hutau_S6=hutau;
clear t huT hutau

load data_S7 t huT hutau
t_S7=t;
huT_S7=huT;
hutau_S7=hutau;
clear t huT hutau





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

norm_huT_S1_sta=huT_S1(Ind_S1:Ind_S1_fi);
norm_huT_S2_sta=huT_S2(Ind_S2:Ind_S2_fi);
norm_huT_S3_sta=huT_S3(Ind_S3:Ind_S3_fi);
norm_huT_S4_sta=huT_S4(Ind_S4:Ind_S4_fi);
norm_huT_S5_sta=huT_S5(Ind_S5:Ind_S5_fi);
norm_huT_S6_sta=huT_S6(Ind_S6:Ind_S6_fi);
norm_huT_S7_sta=huT_S7(Ind_S7:Ind_S7_fi);


MAX_S1=max(abs(norm_huT_S1_sta))
MAX_S2=max(abs(norm_huT_S2_sta))
MAX_S3=max(abs(norm_huT_S3_sta))
MAX_S4=max(abs(norm_huT_S4_sta))
MAX_S5=max(abs(norm_huT_S5_sta))
MAX_S6=max(abs(norm_huT_S6_sta))
MAX_S7=max(abs(norm_huT_S7_sta))

MAX=[MAX_S1;MAX_S2;MAX_S3;MAX_S4;MAX_S5;MAX_S6;MAX_S7]











% %find the index of an event
% indices=find((j-[0;j(1:end-1)])==1);
% %corresponding t value
% t_event=t(indices);
% %inter-event intervals
% t_inter=t_event-[0;t_event(1:end-1)];


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


%===========================================================
axes(h(1))
h7=plot(t_S7,huT_S7,...
    'LineStyle',lines{1},'Color',yellow,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',yellow);hold on;
h6=plot(t_S6,huT_S6,...
    'LineStyle',lines{1},'Color',red,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',red);hold on;
h5=plot(t_S5,huT_S5,...
    'LineStyle',lines{1},'Color',purple,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',purple);hold on;
h4=plot(t_S4,huT_S4,...
    'LineStyle',lines{1},'Color',green,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',green);hold on;
h3=plot(t_S3,huT_S3,...
    'LineStyle',lines{1},'Color',blue,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',blue);hold on;
h2=plot(t_S2,huT_S2,...
    'LineStyle',lines{1},'Color',grey,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',grey);hold on;
h1=plot(t_S1,huT_S1,...
    'LineStyle',lines{1},'Color',black,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',black);hold off;

grid on
xlim([0,100])
% ylim([-320,200])
ylim([-440,440])

leg=legend([h7 h6 h5 h4 h3 h2 h1], {'S7)','S6)','S5)','S4)','S3)','S2)','S1)'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','best')
pos = leg.Position;
set(leg,...
    'Position',pos+[-0.33 -0.3 +0.05 0.35])

leg.Direction = 'reverse';
% leg.ItemTokenSize = [5,5]; 

% uistack(h_aux1, 'top')

% set(gca,'ylim',ylim,'YScale', 'linear')
% yticks([10^(-3) 10^(-2) 10^(-1) 1 10 100]);
% yticklabels({'0.001','0.01','0.1','1','10','100'})
% xlabel('$t\,[s]$','fontsize',12);
% ylabel('$u_T\,[N]$','fontsize',12)
% xh=xlabel('$t\,[s]$','fontsize',12)
yh=ylabel('$u_T\,[N]$','fontsize',12)

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

% pos=get(gca,'position');  % get current axes position vector
% dh=0.05;                   % guess for height adjustment value
% pos(1)=pos(1)+dh;pos(3)=pos(3)-dh;  % raise bottom, reduce height
% set(gca,'position',pos)

%===========================================================
axes(h(2))
h7=plot(t_S7,huT_S7,...
    'LineStyle',lines{1},'Color',yellow,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',yellow);hold on;
h6=plot(t_S6,huT_S6,...
    'LineStyle',lines{1},'Color',red,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',red);hold on;
h5=plot(t_S5,huT_S5,...
    'LineStyle',lines{1},'Color',purple,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',purple);hold on;
h4=plot(t_S4,huT_S4,...
    'LineStyle',lines{1},'Color',green,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',green);hold on;
h3=plot(t_S3,huT_S3,...
    'LineStyle',lines{1},'Color',blue,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',blue);hold on;
h2=plot(t_S2,huT_S2,...
    'LineStyle',lines{1},'Color',grey,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',grey);hold on;
h1=plot(t_S1,huT_S1,...
    'LineStyle',lines{1},'Color',black,'LineWidth',lineWidth,...
    'Marker','none','MarkerSize',markerSize,'MarkerEdgeColor','none','MarkerFaceColor',black);hold off;


grid on
xlim([0,100])
ylim([-160,250])

% leg=legend([h7 h6 h5 h4 h3 h2 h1], {'S7)','S6)','S5)','S4)','S3)','S2)','S1)'},'fontsize',12);
% legend('boxoff')
% leg.Orientation='horizontal';
% set(leg,...
%     'Location','best','NumColumns',7)
% pos = leg.Position;
% set(leg,...
%     'Position',pos+[0.40 -0.11 -0.1 0])

xh=xlabel('$t\,[s]$','fontsize',12)
yh=ylabel('$u_T\,[N]$','fontsize',12)

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

matlabfrag2('simuD_uT')