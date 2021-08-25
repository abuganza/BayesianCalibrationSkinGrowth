clc
clear all;
close all;

%%
Iso_apex=load('Iso_apex_simulation.txt');
Iso_th_apex=load('Iso_apex_th_simulation.txt');
% load theta at day 7
Iso_th_7day=load('Iso_th_7day.txt');
%% Load data
% load 5,50 95 percentile of theta
th_95=Iso_th_apex(1,:);
th_50=Iso_th_apex(2,:);
th_5=Iso_th_apex(3,:);
% load 5,50 95 percentile of theta_g
thg_5=Iso_apex(7,:);
thg_50=Iso_apex(4,:);
thg_95=Iso_apex(1,:);
% calculate 5,50 95 percentile of theta_e
the_95=th_50./thg_5;
the_50=th_50./thg_50;
the_5=th_50./thg_95;


x_fill=[[0:7] fliplr([0:7])];
inBetween_thg =[thg_5, fliplr(thg_95)];
inBetween_th=[th_5, fliplr(th_95)];
inBetween_the=[the_5, fliplr(the_95)];

figure(1)
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis([0 7 1 1.75]);
hold on

fill(x_fill, inBetween_thg, 'r','FaceAlpha',.2,'EdgeColor','none');
fill(x_fill, inBetween_the, 'b','FaceAlpha',.2,'EdgeColor','none');
fill(x_fill, inBetween_th, 'k','FaceAlpha',.2,'EdgeColor','none');
plot([0:7],th_50,'-k','linewidth',2)
plot([0:7],thg_50,'-r','linewidth',2)
plot([0:7],the_50,'-b','linewidth',2)

%%

figure(2)
image(Iso_th_7day,'CDataMapping','scaled');
caxis([1 1.7]);
colormap('jet')
box on
%grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';