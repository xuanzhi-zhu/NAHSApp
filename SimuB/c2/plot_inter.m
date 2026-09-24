clear all
close all

%===============
load data_S1_88 t j zp
t_88=t;
j_88=j;
%find the index of an event
indices=find((j_88-[0;j_88(1:end-1)])==1);
%corresponding t value
t_event_88=t_88(indices);
%inter-event intervals
t_inter_88=t_event_88-[0;t_event_88(1:end-1)];
MIN_88=min(   t_inter_88(2:end));
AVG_88=mean(t_inter_88(2:end));
MAX_88=max( t_inter_88(2:end));
aux=size(zp);
zp_88=reshape(zp,[aux(1) aux(3)]);
norm_zp_88=(sum(zp_88.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_88=max(norm_zp_88(aux:end));
clear t j zp

load data_S1_89 t j zp
t_89=t;
j_89=j;
%find the index of an event
indices=find((j_89-[0;j_89(1:end-1)])==1);
%corresponding t value
t_event_89=t_89(indices);
%inter-event intervals
t_inter_89=t_event_89-[0;t_event_89(1:end-1)];
MIN_89=min(   t_inter_89(2:end));
AVG_89=mean(t_inter_89(2:end));
MAX_89=max( t_inter_89(2:end));
aux=size(zp);
zp_89=reshape(zp,[aux(1) aux(3)]);
norm_zp_89=(sum(zp_89.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_89=max(norm_zp_89(aux:end));
clear t j zp

load data_S1_90 t j zp
t_90=t;
j_90=j;
%find the index of an event
indices=find((j_90-[0;j_90(1:end-1)])==1);
%corresponding t value
t_event_90=t_90(indices);
%inter-event intervals
t_inter_90=t_event_90-[0;t_event_90(1:end-1)];
MIN_90=min(   t_inter_90(2:end));
AVG_90=mean(t_inter_90(2:end));
MAX_90=max( t_inter_90(2:end));
aux=size(zp);
zp_90=reshape(zp,[aux(1) aux(3)]);
norm_zp_90=(sum(zp_90.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_90=max(norm_zp_90(aux:end));
clear t j zp

load data_S1_91 t j zp
t_91=t;
j_91=j;
%find the index of an event
indices=find((j_91-[0;j_91(1:end-1)])==1);
%corresponding t value
t_event_91=t_91(indices);
%inter-event intervals
t_inter_91=t_event_91-[0;t_event_91(1:end-1)];
MIN_91=min(   t_inter_91(2:end));
AVG_91=mean(t_inter_91(2:end));
MAX_91=max( t_inter_91(2:end));
aux=size(zp);
zp_91=reshape(zp,[aux(1) aux(3)]);
norm_zp_91=(sum(zp_91.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_91=max(norm_zp_91(aux:end));
clear t j zp

load data_S1_92 t j zp
t_92=t;
j_92=j;
%find the index of an event
indices=find((j_92-[0;j_92(1:end-1)])==1);
%corresponding t value
t_event_92=t_92(indices);
%inter-event intervals
t_inter_92=t_event_92-[0;t_event_92(1:end-1)];
MIN_92=min(   t_inter_92(2:end));
AVG_92=mean(t_inter_92(2:end));
MAX_92=max( t_inter_92(2:end));
aux=size(zp);
zp_92=reshape(zp,[aux(1) aux(3)]);
norm_zp_92=(sum(zp_92.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_92=max(norm_zp_92(aux:end));
clear t j zp

load data_S1_93 t j zp
t_93=t;
j_93=j;
%find the index of an event
indices=find((j_93-[0;j_93(1:end-1)])==1);
%corresponding t value
t_event_93=t_93(indices);
%inter-event intervals
t_inter_93=t_event_93-[0;t_event_93(1:end-1)];
MIN_93=min(   t_inter_93(2:end));
AVG_93=mean(t_inter_93(2:end));
MAX_93=max( t_inter_93(2:end));
aux=size(zp);
zp_93=reshape(zp,[aux(1) aux(3)]);
norm_zp_93=(sum(zp_93.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_93=max(norm_zp_93(aux:end));
clear t j zp

load data_S1_94 t j zp
t_94=t;
j_94=j;
%find the index of an event
indices=find((j_94-[0;j_94(1:end-1)])==1);
%corresponding t value
t_event_94=t_94(indices);
%inter-event intervals
t_inter_94=t_event_94-[0;t_event_94(1:end-1)];
MIN_94=min(   t_inter_94(2:end));
AVG_94=mean(t_inter_94(2:end));
MAX_94=max( t_inter_94(2:end));
aux=size(zp);
zp_94=reshape(zp,[aux(1) aux(3)]);
norm_zp_94=(sum(zp_94.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_94=max(norm_zp_94(aux:end));
clear t j zp

load data_S1_95 t j zp
t_95=t;
j_95=j;
%f4nd the index of an event
indices=find((j_95-[0;j_95(1:end-1)])==1);
%corresponding t value
t_event_95=t_95(indices);
%inter-event intervals
t_inter_95=t_event_95-[0;t_event_95(1:end-1)];
MIN_95=min(   t_inter_95(2:end));
AVG_95=mean(t_inter_95(2:end));
MAX_95=max( t_inter_95(2:end));
aux=size(zp);
zp_95=reshape(zp,[aux(1) aux(3)]);
norm_zp_95=(sum(zp_95.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_95=max(norm_zp_95(aux:end));
clear t j zp

load data_S1_96 t j zp
t_96=t;
j_96=j;
%find the index of an event
indices=find((j_96-[0;j_96(1:end-1)])==1);
%corresponding t value
t_event_96=t_96(indices);
%inter-event intervals
t_inter_96=t_event_96-[0;t_event_96(1:end-1)];
MIN_96=min(   t_inter_96(2:end));
AVG_96=mean(t_inter_96(2:end));
MAX_96=max( t_inter_96(2:end));
aux=size(zp);
zp_96=reshape(zp,[aux(1) aux(3)]);
norm_zp_96=(sum(zp_96.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_96=max(norm_zp_96(aux:end));
clear t j zp

load data_S1_97 t j zp
t_97=t;
j_97=j;
%find the index of an event
indices=find((j_97-[0;j_97(1:end-1)])==1);
%corresponding t value
t_event_97=t_97(indices);
%inter-event intervals
t_inter_97=t_event_97-[0;t_event_97(1:end-1)];
MIN_97=min(   t_inter_97(2:end));
AVG_97=mean(t_inter_97(2:end));
MAX_97=max( t_inter_97(2:end));
aux=size(zp);
zp_97=reshape(zp,[aux(1) aux(3)]);
norm_zp_97=(sum(zp_97.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_97=max(norm_zp_97(aux:end));
clear t j zp

load data_S1_98 t j zp
t_98=t;
j_98=j;
%find the index of an event
indices=find((j_98-[0;j_98(1:end-1)])==1);
%corresponding t value
t_event_98=t_98(indices);
%inter-event intervals
t_inter_98=t_event_98-[0;t_event_98(1:end-1)];
MIN_98=min(   t_inter_98(2:end));
AVG_98=mean(t_inter_98(2:end));
MAX_98=max( t_inter_98(2:end));
aux=size(zp);
zp_98=reshape(zp,[aux(1) aux(3)]);
norm_zp_98=(sum(zp_98.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_98=max(norm_zp_98(aux:end));
clear t j zp

load data_S1_99 t j zp
t_99=t;
j_99=j;
%find the index of an event
indices=find((j_99-[0;j_99(1:end-1)])==1);
%corresponding t value
t_event_99=t_99(indices);
%inter-event intervals
t_inter_99=t_event_99-[0;t_event_99(1:end-1)];
MIN_99=min(   t_inter_99(2:end));
AVG_99=mean(t_inter_99(2:end));
MAX_99=max( t_inter_99(2:end));
aux=size(zp);
zp_99=reshape(zp,[aux(1) aux(3)]);
norm_zp_99=(sum(zp_99.^2,1)).^(0.5);
aux=min(find(t>=30));
norm_zp_ss_99=max(norm_zp_99(aux:end));
clear t j zp

layout = [1;1;1;1]*ones(1,8);
h=create_axis(layout,18,...
    'innerymargin',0.1,...
    'botmargin',0,...
    'innerxmargin',0.08,...
    'leftmargin',0.02);
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


lineWidth = 1.5;
markerSize = 10;

% aux_k=1.5:1:10.5;

aux_k=0.88:0.01:0.99;

% aux_MIN=[MIN_1_5;MIN_2_5;MIN_3_5;MIN_4_5;MIN_5_5; ...
%                  MIN_6_5;MIN_7_5;MIN_8_5;MIN_9_5;MIN_10_5];
% aux_AVG=[AVG_1_5;AVG_2_5;AVG_3_5;AVG_4_5;AVG_5_5; ...
%                  AVG_6_5;AVG_7_5;AVG_8_5;AVG_9_5;AVG_10_5];
% aux_MAX=[MAX_1_5;MAX_2_5;MAX_3_5;MAX_4_5;MAX_5_5; ...
%                  MAX_6_5;MAX_7_5;MAX_8_5;MAX_9_5;MAX_10_5];

aux_MIN=[MIN_88;MIN_89;MIN_90;MIN_91;MIN_92; ...
                MIN_93;MIN_94;MIN_95;MIN_96;MIN_97;MIN_98;MIN_99];
aux_AVG=[AVG_88;AVG_89;AVG_90;AVG_91;AVG_92; ...
                AVG_93;AVG_94;AVG_95;AVG_96;AVG_97;AVG_98;AVG_99];
aux_MAX=[MAX_88;MAX_89;MAX_90;MAX_91;MAX_92; ...
                MAX_93;MAX_94;MAX_95;MAX_96;MAX_97;MAX_98;MAX_99];

aux_zp_ss=[norm_zp_ss_88;norm_zp_ss_89;...
                   norm_zp_ss_90;norm_zp_ss_91;...
                   norm_zp_ss_92;norm_zp_ss_93;...
                   norm_zp_ss_94;norm_zp_ss_95;...
                   norm_zp_ss_96;norm_zp_ss_97;...
                   norm_zp_ss_98;norm_zp_ss_99];

%===========================================================
axes(h(1))
% colororder([black;grey]);

% yyaxis left
h_auxMIN=plot(aux_k,aux_MIN,...
    markers{8},'MarkerEdgeColor',blue,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_auxAVG=plot(aux_k,aux_AVG,...
    markers{6},'MarkerEdgeColor',black,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_auxMAX=plot(aux_k,aux_MAX,...
    markers{9},'MarkerEdgeColor',red,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_auxss=plot(aux_k,aux_zp_ss,...
    markers{2},'MarkerEdgeColor',grey,'LineWidth',lineWidth,'MarkerSize',markerSize);

set(gca,'YScale', 'log')

% ylabel('$t_{j+1}-t_j\,[s]$','fontsize',12)

grid on

% yyaxis right
% h_auxss=plot(aux_k,aux_zp_ss,...
%     markers{2},'MarkerEdgeColor',grey,'LineWidth',lineWidth,'MarkerSize',markerSize);
% 
% set(gca,'YScale', 'log')
% 
% ylabel('$[m]$','fontsize',12)

grid on
xlim([0.88,0.99])
ylim([0.0001,100])

leg=legend({'$\min_{j\in\naturals}\{t_{j+1}-t_j\}$','$\mathrm{mean}_{j\in\naturals}\{t_{j+1}-t_j\}$','$\max_{j\in\naturals}\{t_{j+1}-t_j\}$','$\max_{t\geq 30}\{|p(t,j)-p_d(t,j)|\}$'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','best','NumColumns',1)
pos = leg.Position;
set(leg,...
    'Position',pos+[-0.03 -0.06 0 0.05])

% ytickLabels = cellstr(num2str(round(log10(yticks(:))), '10^%d'));

% yticks([10^(-4) 10^(-2) 10^(0)]);
% yticklabels({'10^{-4}','10^{-2}','10^{0}'})

xticks([0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99]);
xticklabels({'0.88','0.89','0.90','0.91','0.92','0.93','0.94','0.95','0.96','0.97','0.98','0.99'})

% leg.Direction = 'reverse';
% leg.ItemTokenSize = [5,5]; 

xlabel('$c_2$','fontsize',12);
% ylabel('$t_{j+1}-t_j\,[s]$','fontsize',12)


% 
ax = gca; % current axes
ax.FontSize = 12;
matlabfrag2('simuB_inter_c2')