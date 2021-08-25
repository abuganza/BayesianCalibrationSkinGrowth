clc
clear all;
close all;
%% Read Input
G1_data=readmatrix('G1_Data_Filtered.csv','Range','A2');
G1_data=G1_data';
G1_th=G1_data(1,:);
G1_thg=G1_data(2,:);
G1_thp=G1_data(4,:);
G1_time=G1_data(3,:);

%% Load Pymc3 Result
pymc3_G1=load('pymc3_G1_output_v2.txt');
k_G1=exp(pymc3_G1)/50;
histogram(k_G1,40,'Normalization','probability')
%% G1 Growth vs total deformation 

G1_k_5prc=prctile(k_G1,5);
G1_k_20prc=prctile(k_G1,20);
G1_k_35prc=prctile(k_G1,35);
G1_k_50prc=prctile(k_G1,50);
G1_k_65prc=prctile(k_G1,65);
G1_k_80prc=prctile(k_G1,80);
G1_k_95prc=prctile(k_G1,95);

G1_th_1d=G1_th(1:79);
G1_thg_1d=G1_thg(1:79);
G1_th_3d=G1_th(80:138);
G1_thg_3d=G1_thg(80:138);
G1_th_7d=G1_th(139:end);
G1_thg_7d=G1_thg(139:end);
G1_thg_cal_1d_5=zeros(size(G1_th_1d));
G1_thg_cal_1d_20=zeros(size(G1_th_1d));
G1_thg_cal_1d_35=zeros(size(G1_th_1d));
G1_thg_cal_1d_50=zeros(size(G1_th_1d));
G1_thg_cal_1d_65=zeros(size(G1_th_1d));
G1_thg_cal_1d_80=zeros(size(G1_th_1d));
G1_thg_cal_1d_95=zeros(size(G1_th_1d));
for i=1:length(G1_th_1d)
    temp=find_final_th_g(24,0.1,G1_th_1d(i),G1_thp(1),G1_k_95prc);
    G1_thg_cal_1d_95(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_1d(i),G1_thp(1),G1_k_80prc);
    G1_thg_cal_1d_80(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_1d(i),G1_thp(1),G1_k_65prc);
    G1_thg_cal_1d_65(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_1d(i),G1_thp(1),G1_k_50prc);
    G1_thg_cal_1d_50(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_1d(i),G1_thp(1),G1_k_35prc);
    G1_thg_cal_1d_35(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_1d(i),G1_thp(1),G1_k_20prc);
    G1_thg_cal_1d_20(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_1d(i),G1_thp(1),G1_k_5prc);
    G1_thg_cal_1d_5(i)=temp(end);
end

G1_thg_cal_3d_5=zeros(size(G1_th_3d));
G1_thg_cal_3d_20=zeros(size(G1_th_3d));
G1_thg_cal_3d_35=zeros(size(G1_th_3d));
G1_thg_cal_3d_50=zeros(size(G1_th_3d));
G1_thg_cal_3d_65=zeros(size(G1_th_3d));
G1_thg_cal_3d_80=zeros(size(G1_th_3d));
G1_thg_cal_3d_95=zeros(size(G1_th_3d));
for i=1:length(G1_th_3d)
    temp=find_final_th_g(24,0.1,G1_th_3d(i),G1_thp(1),G1_k_95prc);
    G1_thg_cal_3d_95(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_3d(i),G1_thp(1),G1_k_80prc);
    G1_thg_cal_3d_80(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_3d(i),G1_thp(1),G1_k_65prc);
    G1_thg_cal_3d_65(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_3d(i),G1_thp(1),G1_k_50prc);
    G1_thg_cal_3d_50(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_3d(i),G1_thp(1),G1_k_35prc);
    G1_thg_cal_3d_35(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_3d(i),G1_thp(1),G1_k_20prc);
    G1_thg_cal_3d_20(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_3d(i),G1_thp(1),G1_k_5prc);
    G1_thg_cal_3d_5(i)=temp(end);
end

G1_thg_cal_7d_5=zeros(size(G1_th_7d));
G1_thg_cal_7d_20=zeros(size(G1_th_7d));
G1_thg_cal_7d_35=zeros(size(G1_th_7d));
G1_thg_cal_7d_50=zeros(size(G1_th_7d));
G1_thg_cal_7d_65=zeros(size(G1_th_7d));
G1_thg_cal_7d_80=zeros(size(G1_th_7d));
G1_thg_cal_7d_95=zeros(size(G1_th_7d));
for i=1:length(G1_th_7d)
    temp=find_final_th_g(24,0.1,G1_th_7d(i),G1_thp(1),G1_k_95prc);
    G1_thg_cal_7d_95(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_7d(i),G1_thp(1),G1_k_80prc);
    G1_thg_cal_7d_80(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_7d(i),G1_thp(1),G1_k_65prc);
    G1_thg_cal_7d_65(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_7d(i),G1_thp(1),G1_k_50prc);
    G1_thg_cal_7d_50(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_7d(i),G1_thp(1),G1_k_35prc);
    G1_thg_cal_7d_35(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_7d(i),G1_thp(1),G1_k_20prc);
    G1_thg_cal_7d_20(i)=temp(end);
    temp=find_final_th_g(24,0.1,G1_th_7d(i),G1_thp(1),G1_k_5prc);
    G1_thg_cal_7d_5(i)=temp(end);
end

limit=[1 1.4 1 1.4];
figure(1)
subplot(1,3,1)
hold on
scatter(G1_th_1d,G1_thg_1d,150,'filled','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(G1_th_1d,G1_thg_cal_1d_50,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(G1_th_1d,G1_thg_cal_1d_65,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G1_th_1d,G1_thg_cal_1d_35,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G1_th_1d,G1_thg_cal_1d_80,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G1_th_1d,G1_thg_cal_1d_20,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G1_th_1d,G1_thg_cal_1d_95,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
scatter(G1_th_1d,G1_thg_cal_1d_5,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis(limit)
hold off

subplot(1,3,2)
hold on
scatter(G1_th_3d,G1_thg_3d,150,'filled','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(G1_th_3d,G1_thg_cal_3d_50,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(G1_th_3d,G1_thg_cal_3d_65,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G1_th_3d,G1_thg_cal_3d_35,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G1_th_3d,G1_thg_cal_3d_80,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G1_th_3d,G1_thg_cal_3d_20,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G1_th_3d,G1_thg_cal_3d_95,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
scatter(G1_th_3d,G1_thg_cal_3d_5,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis(limit)
hold off

subplot(1,3,3)
hold on
scatter(G1_th_7d,G1_thg_7d,150,'filled','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(G1_th_7d,G1_thg_cal_7d_50,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(G1_th_7d,G1_thg_cal_7d_65,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G1_th_7d,G1_thg_cal_7d_35,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(G1_th_7d,G1_thg_cal_7d_80,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G1_th_7d,G1_thg_cal_7d_20,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(G1_th_7d,G1_thg_cal_7d_95,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
scatter(G1_th_7d,G1_thg_cal_7d_5,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis(limit)
hold off