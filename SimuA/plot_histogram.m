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



%only steady state inter-transmission time >30 s

T_ss=30;
T_ss_fi=100;
Ind_S1=find(t_event_S1>T_ss, 1 );
Ind_S3=find(t_event_S3>T_ss, 1 );
Ind_S4=find(t_event_S4>T_ss, 1 );
Ind_S5=find(t_event_S5>T_ss, 1 );
Ind_S6=find(t_event_S6>T_ss, 1 );
Ind_S7=find(t_event_S7>T_ss, 1 );

Ind_S1_fi=find(t_event_S1<T_ss_fi, 1,'last' );
Ind_S3_fi=find(t_event_S3<T_ss_fi, 1,'last' );
Ind_S4_fi=find(t_event_S4<T_ss_fi, 1,'last' );
Ind_S5_fi=find(t_event_S5<T_ss_fi, 1,'last' );
Ind_S6_fi=find(t_event_S6<T_ss_fi, 1,'last' );
Ind_S7_fi=find(t_event_S7<T_ss_fi, 1,'last' );


t_event_S1=t_event_S1(Ind_S1:Ind_S1_fi);
t_event_S3=t_event_S3(Ind_S3:Ind_S3_fi);
t_event_S4=t_event_S4(Ind_S4:Ind_S4_fi);
t_event_S5=t_event_S5(Ind_S5:Ind_S5_fi);
t_event_S6=t_event_S6(Ind_S6:Ind_S6_fi);
t_event_S7=t_event_S7(Ind_S7:Ind_S7_fi);

t_inter_S1=t_inter_S1(Ind_S1:Ind_S1_fi);
t_inter_S3=t_inter_S3(Ind_S3:Ind_S3_fi);
t_inter_S4=t_inter_S4(Ind_S4:Ind_S4_fi);
t_inter_S5=t_inter_S5(Ind_S5:Ind_S5_fi);
t_inter_S6=t_inter_S6(Ind_S6:Ind_S6_fi);
t_inter_S7=t_inter_S7(Ind_S7:Ind_S7_fi);



MEAN_S1=mean(t_inter_S1)
MEAN_S3=mean(t_inter_S3)
MEAN_S4=mean(t_inter_S4)
MEAN_S5=mean(t_inter_S5)
MEAN_S6=mean(t_inter_S6)
MEAN_S7=mean(t_inter_S7)

MEAN=[MEAN_S1;MEAN_S3;MEAN_S4;MEAN_S5;MEAN_S6;MEAN_S7]


% layout = [[1;1]*ones(1,3) [2;2]*ones(1,3)];
% h=create_axis(layout,13.4,...
%     'innerymargin',0.05,...
%     'botmargin',0,...
%     'innerxmargin',0.1,...
%     'leftmargin',0.03);

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
green=1/256.*[129 179 101];

purple=1/256.*[151 115 166];

red=1/256.*[185 84 80];
yellow=1/256.*[215 154 2];

grey=[0.1 0.1 0.1];
black=[0 0 0];



lines={'-','--','-.',':'};% lines{1};
markers={'+','o','*','.','x','s','d','^','v','>','<','p','h','|','_'};% markers{1};


lineWidth = 2;
markerSize = 5;

alpha_7 = 0.2; % 0 = fully transparent, 1 = fully opaque
alpha_6 = 0.2;
alpha_5 = 0.2;
alpha_4 = 0.2;
alpha_3 = 0.2;
alpha_1 = 0.2;

%===========================================================

% aux_0=10;
% X1=linspace(-1.5,1,10);
% Y1=linspace(aux_0-6,aux_0-2,10);
% Y2=linspace(aux_0-8,aux_0-4,10);
% t_event_S2=linspace(0,30,10);
% t_inter_S2=aux_0.*ones(numel(t_event_S2),1);



%test
axes(h(1))

h7=histogram(t_inter_S7(2:end),'FaceAlpha',alpha_7);hold on
h6=histogram(t_inter_S6(2:end),'FaceAlpha',alpha_6);hold on
h5=histogram(t_inter_S5(2:end),'FaceAlpha',alpha_5);hold on
h4=histogram(t_inter_S4(2:end),'FaceAlpha',alpha_4);hold on
h3=histogram(t_inter_S3(2:end),'FaceAlpha',alpha_3);hold on
h1=histogram(t_inter_S1(2:end),'FaceAlpha',alpha_1);

% h7.Normalization = 'percentage';
h7.NumBins = 61;
h7.BinWidth = 0.05;
h7.FaceColor = yellow;
h7.EdgeColor = yellow;
% ytickformat("percentage")

% h6.Normalization = 'percentage';
h6.NumBins = 61;
h6.BinWidth = 0.05;
h6.FaceColor = red;
h6.EdgeColor = red;
% ytickformat("percentage")

% h5.Normalization = 'percentage';
h5.NumBins = 61;
h5.BinWidth = 0.05;
h5.FaceColor = purple;
h5.EdgeColor = purple;
% ytickformat("percentage")

% h4.Normalization = 'percentage';
h4.NumBins = 61;
h4.BinWidth = 0.05;
h4.FaceColor = green;
h4.EdgeColor = green;
% ytickformat("percentage")

% h3.Normalization = 'percentage';
h3.NumBins = 61;
h3.BinWidth = 0.05;
h3.FaceColor = blue;
h3.EdgeColor = blue;
% ytickformat("percentage")

% h1.Normalization = 'percentage';
h1.NumBins = 61;
h1.BinWidth = 0.05;
h1.FaceColor = grey;
h1.EdgeColor = grey;
% ytickformat("percentage")




leg=legend([h7 h6 h5 h4 h3 h1], {'S7)','S6)','S5)','S4)','S3)','S1)'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','best')
pos = leg.Position;
set(leg,...
    'Position',pos+[-0.83 -0.3 +0.05 0.35])

leg.Direction = 'reverse';


% % ytickLabels = cellstr(num2str(round(log10(yticks(:))), '10^%d'));
% 
% yticks = ([10^(-4) 10^(-2) 10^(0)]);
% yticklabels({'10^{-4}','10^{-2}','10^{0}'})
% 
% 
% 
% leg.Direction = 'reverse';
% leg.ItemTokenSize = [5,5]; 

% xh=xlabel('$t\,[s]$','fontsize',12)
yh=ylabel('$\text{Counts}$','fontsize',12)

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

ylim([0,100])

% aux_0=10;
% X1=linspace(580-1.5,580+1,10);
% Y1=linspace(aux_0-4,aux_0-1.5,10);
% Y2=linspace(aux_0-6,aux_0-3,10);
% t_event_S2=linspace(580,600,10);
% t_inter_S2=aux_0.*ones(numel(t_event_S2),1);

axes(h(2))

h7=histogram(t_inter_S7(2:end),'FaceAlpha',alpha_7);hold on
h6=histogram(t_inter_S6(2:end),'FaceAlpha',alpha_6);hold on
h5=histogram(t_inter_S5(2:end),'FaceAlpha',alpha_5);hold on
h4=histogram(t_inter_S4(2:end),'FaceAlpha',alpha_4);hold on
h3=histogram(t_inter_S3(2:end),'FaceAlpha',alpha_3);hold on
h1=histogram(t_inter_S1(2:end),'FaceAlpha',alpha_1);hold off

% h7.Normalization = 'percentage';
h7.NumBins = 61;
h7.BinWidth = 0.05;
h7.FaceColor = yellow;
h7.EdgeColor = yellow;
% ytickformat("percentage")

% h6.Normalization = 'percentage';
h6.NumBins = 61;
h6.BinWidth = 0.05;
h6.FaceColor = red;
h6.EdgeColor = red;
% ytickformat("percentage")

% h5.Normalization = 'percentage';
h5.NumBins = 61;
h5.BinWidth = 0.05;
h5.FaceColor = purple;
h5.EdgeColor = purple;
% ytickformat("percentage")

% h4.Normalization = 'percentage';
h4.NumBins = 61;
h4.BinWidth = 0.05;
h4.FaceColor = green;
h4.EdgeColor = green;
% ytickformat("percentage")

% h3.Normalization = 'percentage';
h3.NumBins = 61;
h3.BinWidth = 0.05;
h3.FaceColor = blue;
h3.EdgeColor = blue;
% ytickformat("percentage")

% h1.Normalization = 'percentage';
h1.NumBins = 61;
h1.BinWidth = 0.05;
h1.FaceColor = grey;
h1.EdgeColor = grey;
% ytickformat("percentage")



% 
% leg=legend([h7 h6 h5 h4 h3 h1], {'S7)','S6)','S5)','S4)','S3)','S1)'},'fontsize',12);
% legend('boxoff')
% leg.Orientation='vertical';
% set(leg,...
%     'Location','best')
% pos = leg.Position;
% set(leg,...
%     'Position',pos+[-0.83 -0.3 +0.05 0.35])
% 
% leg.Direction = 'reverse';

% xlim([0,100])
% ylim([-320,200])
ylim([0,20])

xh=xlabel('$t_{j+1}-t_j\,[s]$','fontsize',12)
yh=ylabel('$\text{Counts}$','fontsize',12)

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

matlabfrag2('simuA_inter')