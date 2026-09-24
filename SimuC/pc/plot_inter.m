clear all
close all

T_ss=30;

load data_i_0 t j zp
t_i_0=t;
j_i_0=j;
%find the index of an event
indices=find((j_i_0-[0;j_i_0(1:end-1)])==1);
%corresponding t value
t_event_i_0=t_i_0(indices);
%inter-event intervals
t_inter_i_0=t_event_i_0-[0;t_event_i_0(1:end-1)];
aux=min(find(t_event_i_0>=T_ss));
MIN_i_0=min(   t_inter_i_0(aux:end));
AVG_i_0=mean(t_inter_i_0(aux:end));
MAX_i_0=max( t_inter_i_0(aux:end));
aux=size(zp);
zp_i_0=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_0=(sum(zp_i_0.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_0=max(norm_zp_i_0(aux:end));
clear t j zp

load data_i_1 t j zp
t_i_1=t;
j_i_1=j;
%find the index of an event
indices=find((j_i_1-[0;j_i_1(1:end-1)])==1);
%corresponding t value
t_event_i_1=t_i_1(indices);
%inter-event intervals
t_inter_i_1=t_event_i_1-[0;t_event_i_1(1:end-1)];
aux=min(find(t_event_i_1>=T_ss));
MIN_i_1=min(   t_inter_i_1(aux:end));
AVG_i_1=mean(t_inter_i_1(aux:end));
MAX_i_1=max( t_inter_i_1(aux:end));
aux=size(zp);
zp_i_1=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_1=(sum(zp_i_1.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_1=max(norm_zp_i_1(aux:end));
clear t j zp

load data_i_2 t j zp
t_i_2=t;
j_i_2=j;
%find the index of an event
indices=find((j_i_2-[0;j_i_2(1:end-1)])==1);
%corresponding t value
t_event_i_2=t_i_2(indices);
%inter-event intervals
t_inter_i_2=t_event_i_2-[0;t_event_i_2(1:end-1)];
aux=min(find(t_event_i_2>=T_ss));
MIN_i_2=min(   t_inter_i_2(aux:end));
AVG_i_2=mean(t_inter_i_2(aux:end));
MAX_i_2=max( t_inter_i_2(aux:end));
aux=size(zp);
zp_i_2=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_2=(sum(zp_i_2.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_2=max(norm_zp_i_2(aux:end));
clear t j zp

load data_i_3 t j zp
t_i_3=t;
j_i_3=j;
%find the index of an event
indices=find((j_i_3-[0;j_i_3(1:end-1)])==1);
%corresponding t value
t_event_i_3=t_i_3(indices);
%inter-event intervals
t_inter_i_3=t_event_i_3-[0;t_event_i_3(1:end-1)];
aux=min(find(t_event_i_3>=T_ss));
MIN_i_3=min(   t_inter_i_3(aux:end));
AVG_i_3=mean(t_inter_i_3(aux:end));
MAX_i_3=max( t_inter_i_3(aux:end));
aux=size(zp);
zp_i_3=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_3=(sum(zp_i_3.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_3=max(norm_zp_i_3(aux:end));
clear t j zp

load data_i_4 t j zp
t_i_4=t;
j_i_4=j;
%find the index of an event
indices=find((j_i_4-[0;j_i_4(1:end-1)])==1);
%corresponding t value
t_event_i_4=t_i_4(indices);
%inter-event intervals
t_inter_i_4=t_event_i_4-[0;t_event_i_4(1:end-1)];
aux=min(find(t_event_i_4>=T_ss));
t_inter_i_4=t_inter_i_4(t_inter_i_4>0);
MIN_i_4=min(   t_inter_i_4(aux:end));
AVG_i_4=mean(t_inter_i_4(aux:end));
MAX_i_4=max( t_inter_i_4(aux:end));
aux=size(zp);
zp_i_4=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_4=(sum(zp_i_4.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_4=max(norm_zp_i_4(aux:end));
clear t j zp

load data_i_5 t j zp
t_i_5=t;
j_i_5=j;
%find the index of an event
indices=find((j_i_5-[0;j_i_5(1:end-1)])==1);
%corresponding t value
t_event_i_5=t_i_5(indices);
%inter-event intervals
t_inter_i_5=t_event_i_5-[0;t_event_i_5(1:end-1)];
aux=min(find(t_event_i_5>=T_ss));
t_inter_i_5=t_inter_i_5(t_inter_i_5>0);
MIN_i_5=min(   t_inter_i_5(aux:end));
AVG_i_5=mean(t_inter_i_5(aux:end));
MAX_i_5=max( t_inter_i_5(aux:end));
aux=size(zp);
zp_i_5=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_5=(sum(zp_i_5.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_5=max(norm_zp_i_5(aux:end));
clear t j zp

load data_i_6 t j zp
t_i_6=t;
j_i_6=j;
%find the index of an event
indices=find((j_i_6-[0;j_i_6(1:end-1)])==1);
%corresponding t value
t_event_i_6=t_i_6(indices);
%inter-event intervals
t_inter_i_6=t_event_i_6-[0;t_event_i_6(1:end-1)];
aux=min(find(t_event_i_6>=T_ss));
t_inter_i_6=t_inter_i_6(t_inter_i_6>0);
MIN_i_6=min(   t_inter_i_6(aux:end));
AVG_i_6=mean(t_inter_i_6(aux:end));
MAX_i_6=max( t_inter_i_6(aux:end));
aux=size(zp);
zp_i_6=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_6=(sum(zp_i_6.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_6=max(norm_zp_i_6(aux:end));
clear t j zp

load data_i_7 t j zp
t_i_7=t;
j_i_7=j;
%f4nd the index of an event
indices=find((j_i_7-[0;j_i_7(1:end-1)])==1);
%corresponding t value
t_event_i_7=t_i_7(indices);
%inter-event intervals
t_inter_i_7=t_event_i_7-[0;t_event_i_7(1:end-1)];
aux=min(find(t_event_i_7>=T_ss));
t_inter_i_7=t_inter_i_7(t_inter_i_7>0);
MIN_i_7=min(   t_inter_i_7(aux:end));
AVG_i_7=mean(t_inter_i_7(aux:end));
MAX_i_7=max( t_inter_i_7(aux:end));
aux=size(zp);
zp_i_7=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_7=(sum(zp_i_7.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_7=max(norm_zp_i_7(aux:end));
clear t j zp

load data_i_8 t j zp
t_i_8=t;
j_i_8=j;
%find the index of an event
indices=find((j_i_8-[0;j_i_8(1:end-1)])==1);
%corresponding t value
t_event_i_8=t_i_8(indices);
%inter-event intervals
t_inter_i_8=t_event_i_8-[0;t_event_i_8(1:end-1)];
aux=min(find(t_event_i_8>=T_ss));
t_inter_i_8=t_inter_i_8(t_inter_i_8>0);
MIN_i_8=min(   t_inter_i_8(aux:end));
AVG_i_8=mean(t_inter_i_8(aux:end));
MAX_i_8=max( t_inter_i_8(aux:end));
aux=size(zp);
zp_i_8=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_8=(sum(zp_i_8.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_8=max(norm_zp_i_8(aux:end));
clear t j zp

load data_i_9 t j zp
t_i_9=t;
j_i_9=j;
%find the index of an event
indices=find((j_i_9-[0;j_i_9(1:end-1)])==1);
%corresponding t value
t_event_i_9=t_i_9(indices);
%inter-event intervals
t_inter_i_9=t_event_i_9-[0;t_event_i_9(1:end-1)];
aux=min(find(t_event_i_9>=T_ss));
t_inter_i_9=t_inter_i_9(t_inter_i_9>0);
MIN_i_9=min(   t_inter_i_9(aux:end));
AVG_i_9=mean(t_inter_i_9(aux:end));
MAX_i_9=max( t_inter_i_9(aux:end));
aux=size(zp);
zp_i_9=reshape(zp,[aux(1) aux(3)]);
norm_zp_i_9=(sum(zp_i_9.^2,1)).^(0.5);
aux=min(find(t>=T_ss));
norm_zp_ss_i_9=max(norm_zp_i_9(aux:end));
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

aux_i=0:1:9;

% aux_MIN=[MIN_1_5;MIN_2_5;MIN_3_5;MIN_4_5;MIN_5_5; ...
%                  MIN_6_5;MIN_7_5;MIN_8_5;MIN_9_5;MIN_i_1_5];
% aux_AVG=[AVG_1_5;AVG_2_5;AVG_3_5;AVG_4_5;AVG_5_5; ...
%                  AVG_6_5;AVG_7_5;AVG_8_5;AVG_9_5;AVG_i_1_5];
% aux_MAX=[MAX_1_5;MAX_2_5;MAX_3_5;MAX_4_5;MAX_5_5; ...
%                  MAX_6_5;MAX_7_5;MAX_8_5;MAX_9_5;MAX_i_1_5];

aux_MIN=[MIN_i_0;MIN_i_1;MIN_i_2;MIN_i_3;MIN_i_4; ...
                MIN_i_5;MIN_i_6;MIN_i_7;MIN_i_8;MIN_i_9];
aux_AVG=[AVG_i_0;AVG_i_1;AVG_i_2;AVG_i_3;AVG_i_4; ...
                AVG_i_5;AVG_i_6;AVG_i_7;AVG_i_8;AVG_i_9];
aux_MAX=[MAX_i_0;MAX_i_1;MAX_i_2;MAX_i_3;MAX_i_4; ...
                MAX_i_5;MAX_i_6;MAX_i_7;MAX_i_8;MAX_i_9];

aux_zp_ss=[norm_zp_ss_i_0;norm_zp_ss_i_1;...
                   norm_zp_ss_i_2;norm_zp_ss_i_3;...
                   norm_zp_ss_i_4;norm_zp_ss_i_5;...
                   norm_zp_ss_i_6;norm_zp_ss_i_7;...
                   norm_zp_ss_i_8;norm_zp_ss_i_9];

%===========================================================
axes(h(1))
h_auxMIN=plot(aux_i,aux_MIN,...
    markers{8},'MarkerEdgeColor',blue,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_auxAVG=plot(aux_i,aux_AVG,...
    markers{6},'MarkerEdgeColor',black,'LineWidth',lineWidth,'MarkerSize',markerSize);hold on;
h_auxMAX=plot(aux_i,aux_MAX,...
    markers{9},'MarkerEdgeColor',red,'LineWidth',lineWidth,'MarkerSize',markerSize);hold off;
% h_auxss=plot(aux_i,aux_zp_ss,...
%     markers{2},'MarkerEdgeColor',grey,'LineWidth',lineWidth,'MarkerSize',markerSize);hold off;
grid on

xlim([-1,10])
ylim([1e-6,1])

set(gca,'YScale', 'log')

leg=legend({'$\min_{j\in\naturals}\{t_{j+1}-t_j\}$','$\mathrm{mean}_{j\in\naturals}\{t_{j+1}-t_j\}$'},'fontsize',12);
legend('boxoff')
leg.Orientation='vertical';
set(leg,...
    'Location','best','NumColumns',1)
pos = leg.Position;
set(leg,...
    'Position',pos+[-0.03 -0.01 0 0.05])


% ytickLabels = cellstr(num2str(round(log10(yticks(:))), '10^%d'));

% yticks = ([10^(-4) 10^(-2) 10^(0)]);
% yticklabels({'10^{-4}','10^{-2}','10^{0}'})



% leg.Direction = 'reverse';
% leg.ItemTokenSize = [5,5]; 

xlabel('$i$','fontsize',12);
ylabel('$t_{j+1}-t_j\,[s]$','fontsize',12)



ax = gca; % current axes
ax.FontSize = 12;
matlabfrag2('simuC_inter_i')