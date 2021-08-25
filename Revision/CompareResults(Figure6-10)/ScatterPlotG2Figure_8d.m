clc
clear all;
close all;
%% Read Input
G2_data=readmatrix('G2_Data_Filtered.csv','Range','A2');
G2_data=G2_data';
G2_th=G2_data(1,:);
G2_thg=G2_data(2,:);
G2_thp=G2_data(4,:);
G2_time=G2_data(3,:);

%% Load Pymc3 Result
pymc3_G2=load('pymc3_G2_output_v2.txt');
k_G2=exp(pymc3_G2)/0.5;
histogram(k_G2,40,'Normalization','probability')
%% G2 Growth vs total deformation 

G2_k_5prc=prctile(k_G2,5);
G2_k_20prc=prctile(k_G2,20);
G2_k_35prc=prctile(k_G2,35);
G2_k_50prc=prctile(k_G2,50);
G2_k_65prc=prctile(k_G2,65);
G2_k_80prc=prctile(k_G2,80);
G2_k_95prc=prctile(k_G2,95);

G2_th_1d=G2_th(1:37);
G2_thg_1d=G2_thg(1:37);
G2_th_3d=G2_th(38:69);
G2_thg_3d=G2_thg(38:69);
G2_th_7d=G2_th(70:end);
G2_thg_7d=G2_thg(70:end);
G2_thg_cal_1d_5=zeros(size(G2_th_1d));
G2_thg_cal_1d_20=zeros(size(G2_th_1d));
G2_thg_cal_1d_35=zeros(size(G2_th_1d));
G2_thg_cal_1d_50=zeros(size(G2_th_1d));
G2_thg_cal_1d_65=zeros(size(G2_th_1d));
G2_thg_cal_1d_80=zeros(size(G2_th_1d));
G2_thg_cal_1d_95=zeros(size(G2_th_1d));
for i=1:length(G2_th_1d)
    temp=find_final_th_g(24,0.1,G2_th_1d(i),G2_thp(1),G2_k_95prc);
    G2_thg_cal_1d_95(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_1d(i),G2_thp(1),G2_k_80prc);
    G2_thg_cal_1d_80(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_1d(i),G2_thp(1),G2_k_65prc);
    G2_thg_cal_1d_65(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_1d(i),G2_thp(1),G2_k_50prc);
    G2_thg_cal_1d_50(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_1d(i),G2_thp(1),G2_k_35prc);
    G2_thg_cal_1d_35(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_1d(i),G2_thp(1),G2_k_20prc);
    G2_thg_cal_1d_20(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_1d(i),G2_thp(1),G2_k_5prc);
    G2_thg_cal_1d_5(i)=temp(end);
end

G2_thg_cal_3d_5=zeros(size(G2_th_3d));
G2_thg_cal_3d_20=zeros(size(G2_th_3d));
G2_thg_cal_3d_35=zeros(size(G2_th_3d));
G2_thg_cal_3d_50=zeros(size(G2_th_3d));
G2_thg_cal_3d_65=zeros(size(G2_th_3d));
G2_thg_cal_3d_80=zeros(size(G2_th_3d));
G2_thg_cal_3d_95=zeros(size(G2_th_3d));
for i=1:length(G2_th_3d)
    temp=find_final_th_g(24,0.1,G2_th_3d(i),G2_thp(1),G2_k_95prc);
    G2_thg_cal_3d_95(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_3d(i),G2_thp(1),G2_k_80prc);
    G2_thg_cal_3d_80(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_3d(i),G2_thp(1),G2_k_65prc);
    G2_thg_cal_3d_65(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_3d(i),G2_thp(1),G2_k_50prc);
    G2_thg_cal_3d_50(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_3d(i),G2_thp(1),G2_k_35prc);
    G2_thg_cal_3d_35(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_3d(i),G2_thp(1),G2_k_20prc);
    G2_thg_cal_3d_20(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_3d(i),G2_thp(1),G2_k_5prc);
    G2_thg_cal_3d_5(i)=temp(end);
end

G2_thg_cal_7d_5=zeros(size(G2_th_7d));
G2_thg_cal_7d_20=zeros(size(G2_th_7d));
G2_thg_cal_7d_35=zeros(size(G2_th_7d));
G2_thg_cal_7d_50=zeros(size(G2_th_7d));
G2_thg_cal_7d_65=zeros(size(G2_th_7d));
G2_thg_cal_7d_80=zeros(size(G2_th_7d));
G2_thg_cal_7d_95=zeros(size(G2_th_7d));
for i=1:length(G2_th_7d)
    temp=find_final_th_g(24,0.1,G2_th_7d(i),G2_thp(1),G2_k_95prc);
    G2_thg_cal_7d_95(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_7d(i),G2_thp(1),G2_k_80prc);
    G2_thg_cal_7d_80(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_7d(i),G2_thp(1),G2_k_65prc);
    G2_thg_cal_7d_65(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_7d(i),G2_thp(1),G2_k_50prc);
    G2_thg_cal_7d_50(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_7d(i),G2_thp(1),G2_k_35prc);
    G2_thg_cal_7d_35(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_7d(i),G2_thp(1),G2_k_20prc);
    G2_thg_cal_7d_20(i)=temp(end);
    temp=find_final_th_g(24,0.1,G2_th_7d(i),G2_thp(1),G2_k_5prc);
    G2_thg_cal_7d_5(i)=temp(end);
end

limit=[1 1.5 1 1.5];
figure(1)
subplot(1,3,1)
hold on
scatter(G2_th_1d,G2_thg_1d,150,'filled','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(G2_th_1d,G2_thg_cal_1d_50,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(G2_th_1d,G2_thg_cal_1d_65,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G2_th_1d,G2_thg_cal_1d_35,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G2_th_1d,G2_thg_cal_1d_80,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G2_th_1d,G2_thg_cal_1d_20,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G2_th_1d,G2_thg_cal_1d_95,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
scatter(G2_th_1d,G2_thg_cal_1d_5,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis(limit);
hold off

subplot(1,3,2)
hold on
scatter(G2_th_3d,G2_thg_3d,150,'filled','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(G2_th_3d,G2_thg_cal_3d_50,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(G2_th_3d,G2_thg_cal_3d_65,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G2_th_3d,G2_thg_cal_3d_35,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G2_th_3d,G2_thg_cal_3d_80,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G2_th_3d,G2_thg_cal_3d_20,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G2_th_3d,G2_thg_cal_3d_95,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
scatter(G2_th_3d,G2_thg_cal_3d_5,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis(limit)
hold off

subplot(1,3,3)
hold on
scatter(G2_th_7d,G2_thg_7d,150,'filled','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(G2_th_7d,G2_thg_cal_7d_50,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(G2_th_7d,G2_thg_cal_7d_65,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G2_th_7d,G2_thg_cal_7d_35,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G2_th_7d,G2_thg_cal_7d_80,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G2_th_7d,G2_thg_cal_7d_20,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G2_th_7d,G2_thg_cal_7d_95,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
scatter(G2_th_7d,G2_thg_cal_7d_5,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis(limit)
hold off