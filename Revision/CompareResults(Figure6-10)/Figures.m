clc
clear all;
close all;

%% Read Input
Iso_data=readmatrix('Iso_Data_Filtered.csv','Range','A2');
G1_data=readmatrix('G1_Data_Filtered.csv','Range','A2');
G2_data=readmatrix('G2_Data_Filtered.csv','Range','A2');
Iso_data=Iso_data';
G1_data=G1_data';
G2_data=G2_data';

%% Read th, thg, thp, time in Iso, G1 and G2
Iso_th=Iso_data(1,:);
Iso_thg=Iso_data(2,:);
Iso_thp=Iso_data(4,:);
Iso_time=Iso_data(3,:);

G1_th=G1_data(1,:);
G1_thg=G1_data(2,:);
G1_thp=G1_data(4,:);
G1_time=G1_data(3,:);

G2_th=G2_data(1,:);
G2_thg=G2_data(2,:);
G2_thp=G2_data(4,:);
G2_time=G2_data(3,:);

%% Process Pymc3 ouput to get growth parameter
% load output
pymc3_Iso=load('pymc3_Iso_output.txt');
%pymc3_G1=load('pymc3_G1_output.txt');
%pymc3_G2=load('pymc3_G2_output.txt');
% calculate growth parameter
k_iso=exp(pymc3_Iso)/25;
%k_G1=exp(pymc3_G1)/50;
%k_G2=exp(pymc3_G2)/0.5;
% %% Histogram of growth parameters
% figure(4)
% subplot(1,3,1)
% hold on
% histogram(k_iso,40,'Normalization','Probability','EdgeColor','k','FaceColor','k','FaceAlpha',0.4)
% xlabel('k_{Iso}');
% ylabel('Probability');
% ax=gca;
% ax.FontSize = 12;
% %set(gcf, 'Position',  [100, 100, 600, 600])
% box on
% hold off
% subplot(1,3,2)
% hold on
% histogram(k_G1,40,'Normalization','Probability','EdgeColor','k','FaceColor','k','FaceAlpha',0.4)
% xlabel('k_{G_1}');
% ylabel('Probability');
% ax=gca;
% ax.FontSize = 12;
% %set(gcf, 'Position',  [100, 100, 600, 600])
% box on
% hold off
% subplot(1,3,3)
% hold on
% histogram(k_G2,40,'Normalization','Probability','EdgeColor','k','FaceColor','k','FaceAlpha',0.4)
% xlabel('k_{G_2}');
% ylabel('Probability');
% ax=gca;
% ax.FontSize = 12;
% %set(gcf, 'Position',  [100, 100, 600, 600])
% box on
% hold off

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
figure(5)
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
hold off
%%
% k_5=0.04;
% k_95=0.041;
% ap_1=10;
% ap_3=94;
% ap_7=149;
% [th_g_cla_ap1_lo]=find_final_th_g(7*24,1,Iso_th(ap_1),Iso_thp(ap_1),0.03);
% [th_g_cla_ap1_up]=find_final_th_g(7*24,1,Iso_th(ap_1),Iso_thp(ap_1),0.05);
% [th_g_cla_ap3_lo]=find_final_th_g(7*24,1,Iso_th(ap_3),Iso_thp(ap_3),0.03);
% [th_g_cla_ap3_up]=find_final_th_g(7*24,1,Iso_th(ap_3),Iso_thp(ap_3),0.05);
% [th_g_cla_ap7_lo]=find_final_th_g(7*24,1,Iso_th(ap_7),Iso_thp(ap_7),0.03);
% [th_g_cla_ap7_up]=find_final_th_g(7*24,1,Iso_th(ap_7),Iso_thp(ap_7),0.05);
% 
% figure(5)
% hold on
% plot([0:1:7*24]/24,th_g_cla_ap1_lo);
% plot([0:1:7*24]/24,th_g_cla_ap1_up);
% plot([0:1:7*24]/24,th_g_cla_ap3_lo);
% plot([0:1:7*24]/24,th_g_cla_ap3_up);
% plot([0:1:7*24]/24,th_g_cla_ap7_lo);
% plot([0:1:7*24]/24,th_g_cla_ap7_up);
% scatter([1,3,7],[Iso_thg(ap_1),Iso_thg(ap_3),Iso_thg(ap_7)],'k','filled');
% legend('ap1_lo','ap1_up','ap3_lo','ap3_up','ap7_lo','ap7_up');