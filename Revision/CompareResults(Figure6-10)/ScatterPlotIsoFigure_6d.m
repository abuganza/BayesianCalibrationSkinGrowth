clc
clear all;
close all;
%% Read Input
Iso_data=readmatrix('Iso_Data_Filtered.csv','Range','A2');
Iso_data=Iso_data';
Iso_th=Iso_data(1,:);
Iso_thg=Iso_data(2,:);
Iso_thp=Iso_data(4,:);
Iso_time=Iso_data(3,:);

%% Load Pymc3 Result
pymc3_Iso=load('pymc3_Iso_output_v2.txt');
k_iso=exp(pymc3_Iso)/25;
histogram(k_iso,40,'Normalization','probability')
%% Iso Growth vs total deformation 
Iso_k_5prc=prctile(k_iso,5);
Iso_k_20prc=prctile(k_iso,20);
Iso_k_35prc=prctile(k_iso,35);
Iso_k_50prc=prctile(k_iso,50);
Iso_k_65prc=prctile(k_iso,65);
Iso_k_80prc=prctile(k_iso,80);
Iso_k_95prc=prctile(k_iso,95);

Iso_th_1d=Iso_th(1:67);
Iso_thg_1d=Iso_thg(1:67);
Iso_th_3d=Iso_th(68:112);
Iso_thg_3d=Iso_thg(68:112);
Iso_th_7d=Iso_th(113:end);
Iso_thg_7d=Iso_thg(113:end);
Iso_thg_cal_1d_5=zeros(size(Iso_th_1d));
Iso_thg_cal_1d_20=zeros(size(Iso_th_1d));
Iso_thg_cal_1d_35=zeros(size(Iso_th_1d));
Iso_thg_cal_1d_50=zeros(size(Iso_th_1d));
Iso_thg_cal_1d_65=zeros(size(Iso_th_1d));
Iso_thg_cal_1d_80=zeros(size(Iso_th_1d));
Iso_thg_cal_1d_95=zeros(size(Iso_th_1d));
for i=1:length(Iso_th_1d)
    temp=find_final_th_g(24,0.1,Iso_th_1d(i),Iso_thp(1),Iso_k_95prc);
    Iso_thg_cal_1d_95(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_1d(i),Iso_thp(1),Iso_k_80prc);
    Iso_thg_cal_1d_80(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_1d(i),Iso_thp(1),Iso_k_65prc);
    Iso_thg_cal_1d_65(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_1d(i),Iso_thp(1),Iso_k_50prc);
    Iso_thg_cal_1d_50(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_1d(i),Iso_thp(1),Iso_k_35prc);
    Iso_thg_cal_1d_35(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_1d(i),Iso_thp(1),Iso_k_20prc);
    Iso_thg_cal_1d_20(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_1d(i),Iso_thp(1),Iso_k_5prc);
    Iso_thg_cal_1d_5(i)=temp(end);
end

Iso_thg_cal_3d_5=zeros(size(Iso_th_3d));
Iso_thg_cal_3d_20=zeros(size(Iso_th_3d));
Iso_thg_cal_3d_35=zeros(size(Iso_th_3d));
Iso_thg_cal_3d_50=zeros(size(Iso_th_3d));
Iso_thg_cal_3d_65=zeros(size(Iso_th_3d));
Iso_thg_cal_3d_80=zeros(size(Iso_th_3d));
Iso_thg_cal_3d_95=zeros(size(Iso_th_3d));
for i=1:length(Iso_th_3d)
    temp=find_final_th_g(24,0.1,Iso_th_3d(i),Iso_thp(1),Iso_k_95prc);
    Iso_thg_cal_3d_95(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_3d(i),Iso_thp(1),Iso_k_80prc);
    Iso_thg_cal_3d_80(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_3d(i),Iso_thp(1),Iso_k_65prc);
    Iso_thg_cal_3d_65(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_3d(i),Iso_thp(1),Iso_k_50prc);
    Iso_thg_cal_3d_50(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_3d(i),Iso_thp(1),Iso_k_35prc);
    Iso_thg_cal_3d_35(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_3d(i),Iso_thp(1),Iso_k_20prc);
    Iso_thg_cal_3d_20(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_3d(i),Iso_thp(1),Iso_k_5prc);
    Iso_thg_cal_3d_5(i)=temp(end);
end

Iso_thg_cal_7d_5=zeros(size(Iso_th_7d));
Iso_thg_cal_7d_20=zeros(size(Iso_th_7d));
Iso_thg_cal_7d_35=zeros(size(Iso_th_7d));
Iso_thg_cal_7d_50=zeros(size(Iso_th_7d));
Iso_thg_cal_7d_65=zeros(size(Iso_th_7d));
Iso_thg_cal_7d_80=zeros(size(Iso_th_7d));
Iso_thg_cal_7d_95=zeros(size(Iso_th_7d));
for i=1:length(Iso_th_7d)
    temp=find_final_th_g(24,0.1,Iso_th_7d(i),Iso_thp(1),Iso_k_95prc);
    Iso_thg_cal_7d_95(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_7d(i),Iso_thp(1),Iso_k_80prc);
    Iso_thg_cal_7d_80(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_7d(i),Iso_thp(1),Iso_k_65prc);
    Iso_thg_cal_7d_65(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_7d(i),Iso_thp(1),Iso_k_50prc);
    Iso_thg_cal_7d_50(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_7d(i),Iso_thp(1),Iso_k_35prc);
    Iso_thg_cal_7d_35(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_7d(i),Iso_thp(1),Iso_k_20prc);
    Iso_thg_cal_7d_20(i)=temp(end);
    temp=find_final_th_g(24,0.1,Iso_th_7d(i),Iso_thp(1),Iso_k_5prc);
    Iso_thg_cal_7d_5(i)=temp(end);
end

limit=[1 1.6 1 1.6];

figure(1)
subplot(1,3,1)
hold on
scatter(Iso_th_1d,Iso_thg_1d,150,'filled','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(Iso_th_1d,Iso_thg_cal_1d_50,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(Iso_th_1d,Iso_thg_cal_1d_65,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(Iso_th_1d,Iso_thg_cal_1d_35,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(Iso_th_1d,Iso_thg_cal_1d_80,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(Iso_th_1d,Iso_thg_cal_1d_20,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(Iso_th_1d,Iso_thg_cal_1d_95,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
scatter(Iso_th_1d,Iso_thg_cal_1d_5,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis(limit);
hold off

subplot(1,3,2)
hold on
scatter(Iso_th_3d,Iso_thg_3d,150,'filled','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(Iso_th_3d,Iso_thg_cal_3d_50,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(Iso_th_3d,Iso_thg_cal_3d_65,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(Iso_th_3d,Iso_thg_cal_3d_35,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(Iso_th_3d,Iso_thg_cal_3d_80,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(Iso_th_3d,Iso_thg_cal_3d_20,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(Iso_th_3d,Iso_thg_cal_3d_95,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
scatter(Iso_th_3d,Iso_thg_cal_3d_5,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis(limit)
hold off

subplot(1,3,3)
hold on
scatter(Iso_th_7d,Iso_thg_7d,150,'filled','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(Iso_th_7d,Iso_thg_cal_7d_50,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(Iso_th_7d,Iso_thg_cal_7d_65,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(Iso_th_7d,Iso_thg_cal_7d_35,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.4);
scatter(Iso_th_7d,Iso_thg_cal_7d_80,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(Iso_th_7d,Iso_thg_cal_7d_20,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(Iso_th_7d,Iso_thg_cal_7d_95,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
scatter(Iso_th_7d,Iso_thg_cal_7d_5,150,'filled','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis(limit)
hold off