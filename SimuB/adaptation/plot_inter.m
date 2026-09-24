clear all
close all

T_ss=30;

load data_scale_m_4 t j zp
t_m_4=t;
j_m_4=j;
%find the index of an event
indices=find((j_m_4-[0;j_m_4(1:end-1)])==1);
%corresponding t value
t_event_m_4=t_m_4(indices);
%inter-event intervals
t_inter_m_4=t_event_m_4-[0;t_event_m_4(1:end-1)];
aux=min(find(t_event_m_4>=T_ss));
MIN_m_4=min(   t_inter_m_4(aux:end));
AVG_m_4=mean(t_inter_m_4(aux:end));
MAX_m_4=max( t_inter_m_4(aux:end));
aux=size(zp);
zp_m_4=reshape(zp,[aux(1) aux(3)]);
norm_zp_m_4=(sum(zp_m_4.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_m_4=max(norm_zp_m_4(aux:end));
clear t j zp

load data_scale_m_3 t j zp
t_m_3=t;
j_m_3=j;
%find the index of an event
indices=find((j_m_3-[0;j_m_3(1:end-1)])==1);
%corresponding t value
t_event_m_3=t_m_3(indices);
%inter-event intervals
t_inter_m_3=t_event_m_3-[0;t_event_m_3(1:end-1)];
aux=min(find(t_event_m_3>=T_ss));
MIN_m_3=min(   t_inter_m_3(aux:end));
AVG_m_3=mean(t_inter_m_3(aux:end));
MAX_m_3=max( t_inter_m_3(aux:end));
aux=size(zp);
zp_m_3=reshape(zp,[aux(1) aux(3)]);
norm_zp_m_3=(sum(zp_m_3.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_m_3=max(norm_zp_m_3(aux:end));
clear t j zp

load data_scale_m_2 t j zp
t_m_2=t;
j_m_2=j;
%find the index of an event
indices=find((j_m_2-[0;j_m_2(1:end-1)])==1);
%corresponding t value
t_event_m_2=t_m_2(indices);
%inter-event intervals
t_inter_m_2=t_event_m_2-[0;t_event_m_2(1:end-1)];
aux=min(find(t_event_m_2>=T_ss));
MIN_m_2=min(   t_inter_m_2(aux:end));
AVG_m_2=mean(t_inter_m_2(aux:end));
MAX_m_2=max( t_inter_m_2(aux:end));
aux=size(zp);
zp_m_2=reshape(zp,[aux(1) aux(3)]);
norm_zp_m_2=(sum(zp_m_2.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_m_2=max(norm_zp_m_2(aux:end));
clear t j zp

load data_scale_m_1 t j zp
t_m_1=t;
j_m_1=j;
%find the index of an event
indices=find((j_m_1-[0;j_m_1(1:end-1)])==1);
%corresponding t value
t_event_m_1=t_m_1(indices);
%inter-event intervals
t_inter_m_1=t_event_m_1-[0;t_event_m_1(1:end-1)];
aux=min(find(t_event_m_1>=T_ss));
MIN_m_1=min(   t_inter_m_1(aux:end));
AVG_m_1=mean(t_inter_m_1(aux:end));
MAX_m_1=max( t_inter_m_1(aux:end));
aux=size(zp);
zp_m_1=reshape(zp,[aux(1) aux(3)]);
norm_zp_m_1=(sum(zp_m_1.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_m_1=max(norm_zp_m_1(aux:end));
clear t j zp

load data_scale_p_0 t j zp
t_p_0=t;
j_p_0=j;
%find the index of an event
indices=find((j_p_0-[0;j_p_0(1:end-1)])==1);
%corresponding t value
t_event_p_0=t_p_0(indices);
%inter-event intervals
t_inter_p_0=t_event_p_0-[0;t_event_p_0(1:end-1)];
aux=min(find(t_event_p_0>=T_ss));
MIN_p_0=min(   t_inter_p_0(aux:end));
AVG_p_0=mean(t_inter_p_0(aux:end));
MAX_p_0=max( t_inter_p_0(aux:end));
aux=size(zp);
zp_p_0=reshape(zp,[aux(1) aux(3)]);
norm_zp_p_0=(sum(zp_p_0.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_p_0=max(norm_zp_p_0(aux:end));
clear t j zp

load data_scale_p_1 t j zp
t_p_1=t;
j_p_1=j;
%find the index of an event
indices=find((j_p_1-[0;j_p_1(1:end-1)])==1);
%corresponding t value
t_event_p_1=t_p_1(indices);
%inter-event intervals
t_inter_p_1=t_event_p_1-[0;t_event_p_1(1:end-1)];
aux=min(find(t_event_p_1>=T_ss));
MIN_p_1=min(   t_inter_p_1(aux:end));
AVG_p_1=mean(t_inter_p_1(aux:end));
MAX_p_1=max( t_inter_p_1(aux:end));
aux=size(zp);
zp_p_1=reshape(zp,[aux(1) aux(3)]);
norm_zp_p_1=(sum(zp_p_1.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_p_1=max(norm_zp_p_1(aux:end));
clear t j zp

load data_scale_p_2 t j zp
t_p_2=t;
j_p_2=j;
%f4nd the index of an event
indices=find((j_p_2-[0;j_p_2(1:end-1)])==1);
%corresponding t value
t_event_p_2=t_p_2(indices);
%inter-event intervals
t_inter_p_2=t_event_p_2-[0;t_event_p_2(1:end-1)];
aux=min(find(t_event_p_2>=T_ss));
MIN_p_2=min(   t_inter_p_2(aux:end));
AVG_p_2=mean(t_inter_p_2(aux:end));
MAX_p_2=max( t_inter_p_2(aux:end));
aux=size(zp);
zp_p_2=reshape(zp,[aux(1) aux(3)]);
norm_zp_p_2=(sum(zp_p_2.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_p_2=max(norm_zp_p_2(aux:end));
clear t j zp

load data_scale_p_3 t j zp
t_p_3=t;
j_p_3=j;
%find the index of an event
indices=find((j_p_3-[0;j_p_3(1:end-1)])==1);
%corresponding t value
t_event_p_3=t_p_3(indices);
%inter-event intervals
t_inter_p_3=t_event_p_3-[0;t_event_p_3(1:end-1)];
aux=min(find(t_event_p_3>=T_ss));
MIN_p_3=min(   t_inter_p_3(aux:end));
AVG_p_3=mean(t_inter_p_3(aux:end));
MAX_p_3=max( t_inter_p_3(aux:end));
aux=size(zp);
zp_p_3=reshape(zp,[aux(1) aux(3)]);
norm_zp_p_3=(sum(zp_p_3.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_p_3=max(norm_zp_p_3(aux:end));
clear t j zp

load data_scale_p_4 t j zp
t_p_4=t;
j_p_4=j;
%find the index of an event
indices=find((j_p_4-[0;j_p_4(1:end-1)])==1);
%corresponding t value
t_event_p_4=t_p_4(indices);
%inter-event intervals
t_inter_p_4=t_event_p_4-[0;t_event_p_4(1:end-1)];
aux=min(find(t_event_p_4>=T_ss));
MIN_p_4=min(   t_inter_p_4(aux:end));
AVG_p_4=mean(t_inter_p_4(aux:end));
MAX_p_4=max( t_inter_p_4(aux:end));
aux=size(zp);
zp_p_4=reshape(zp,[aux(1) aux(3)]);
norm_zp_p_4=(sum(zp_p_4.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_p_4=max(norm_zp_p_4(aux:end));
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

aux_k=-4:1:4;

% aux_MIN=[MIN_1_5;MIN_2_5;MIN_3_5;MIN_4_5;MIN_5_5; ...
%                  MIN_6_5;MIN_7_5;MIN_8_5;MIN_9_5;MIN_m_4_5];
% aux_AVG=[AVG_1_5;AVG_2_5;AVG_3_5;AVG_4_5;AVG_5_5; ...
%                  AVG_6_5;AVG_7_5;AVG_8_5;AVG_9_5;AVG_m_4_5];
% aux_MAX=[MAX_1_5;MAX_2_5;MAX_3_5;MAX_4_5;MAX_5_5; ...
%                  MAX_6_5;MAX_7_5;MAX_8_5;MAX_9_5;MAX_m_4_5];

aux_MIN=[MIN_m_4;MIN_m_3;MIN_m_2;MIN_m_1; ...
                MIN_p_0;MIN_p_1;MIN_p_2;MIN_p_3;MIN_p_4];
aux_AVG=[AVG_m_4;AVG_m_3;AVG_m_2;AVG_m_1; ...
                AVG_p_0;AVG_p_1;AVG_p_2;AVG_p_3;AVG_p_4];
aux_MAX=[MAX_m_4;MAX_m_3;MAX_m_2;MAX_m_1; ...
                MAX_p_0;MAX_p_1;MAX_p_2;MAX_p_3;MAX_p_4];

aux_zp_ss=[norm_zp_ss_m_4;...
                   norm_zp_ss_m_3;norm_zp_ss_m_2;...
                   norm_zp_ss_m_1;norm_zp_ss_p_0;...
                   norm_zp_ss_p_1;norm_zp_ss_p_2;...
                   norm_zp_ss_p_3;norm_zp_ss_p_4];

%===========================================================
axes(h(1))
h_auxMIN=plot(aux_k,aux_MIN,...
    markers{8},'MarkerEdgeColor',blue,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_auxAVG=plot(aux_k,aux_AVG,...
    markers{6},'MarkerEdgeColor',black,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_auxMAX=plot(aux_k,aux_MAX,...
    markers{9},'MarkerEdgeColor',red,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_auxss=plot(aux_k,aux_zp_ss,...
    markers{2},'MarkerEdgeColor',grey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold off;
grid on

xlim([-5,5])
ylim([0.1,19])

set(gca,'YScale', 'log')

leg=legend({'$\min_{j\in\naturals}\{t_{j+1}-t_j\}$','$\mathrm{mean}_{j\in\naturals}\{t_{j+1}-t_j\}$','$\max_{j\in\naturals}\{t_{j+1}-t_j\}$','$\max_{t\geq 30}\{|p(t,j)-p_d(t,j)|\}$'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','best','NumColumns',1)
pos = leg.Position;
set(leg,...
    'Position',pos+[-0.03 -0.08 0 0.05])


% ytickLabels = cellstr(num2str(round(log10(yticks(:))), '10^%d'));

% yticks = ([10^(-4) 10^(-2) 10^(0)]);
% yticklabels({'10^{-4}','10^{-2}','10^{0}'})



% leg.Direction = 'reverse';
% leg.ItemTokenSize = [5,5]; 

xlabel('$k$','fontsize',12);
% ylabel('$t_{j+1}-t_j\,[s]$','fontsize',12)



ax = gca; % current axes
ax.FontSize = 12;
matlabfrag2('simuB_inter_adaptation')