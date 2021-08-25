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
histogram(k_iso,40,'Normalization','probability','FaceColor',[0.4 0.4 0.4]);
Iso_k_5prc=prctile(k_iso,5);
Iso_k_20prc=prctile(k_iso,20);
Iso_k_35prc=prctile(k_iso,35);
Iso_k_50prc=prctile(k_iso,50);
Iso_k_65prc=prctile(k_iso,65);
Iso_k_80prc=prctile(k_iso,80);
Iso_k_95prc=prctile(k_iso,95);

%% plot of Iso th_g over time
ap_1=15;
ap_3=81;
ap_7=152;
figure(1)
hold on
box on
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
hold off

%% Plot of th_g over time with DE siumlation

% index of the 9 data points used
ap_1=41;
ap_3=102;
ap_7=147;

in_1=8;
in_3=72; 
in_7=285;

pe_1=3;
pe_3=78;
pe_7=163;

% calculate upper and lower bound of apex
[th_g_ap_95]=find_final_th_g(7*24,1,Iso_th(ap_7),Iso_thp(ap_7),Iso_k_95prc);
[th_g_ap_80]=find_final_th_g(7*24,1,Iso_th(ap_7),Iso_thp(ap_7),Iso_k_80prc);
[th_g_ap_65]=find_final_th_g(7*24,1,Iso_th(ap_7),Iso_thp(ap_7),Iso_k_65prc);
[th_g_ap_5]=find_final_th_g(7*24,1,Iso_th(ap_1),Iso_thp(ap_1),Iso_k_5prc);
[th_g_ap_20]=find_final_th_g(7*24,1,Iso_th(ap_1),Iso_thp(ap_1),Iso_k_20prc);
[th_g_ap_35]=find_final_th_g(7*24,1,Iso_th(ap_1),Iso_thp(ap_1),Iso_k_35prc);

% mark out shaded area for apex
x_fill=[[0:7*24]/24 fliplr([0:7*24]/24)];
inBetween_apex_90=[th_g_ap_5, fliplr(th_g_ap_95)];
inBetween_apex_60=[th_g_ap_20, fliplr(th_g_ap_80)];
inBetween_apex_30=[th_g_ap_35, fliplr(th_g_ap_65)];

% calculate upper and lower bound of intermidiate
[th_g_in_95]=find_final_th_g(7*24,1,Iso_th(in_7),Iso_thp(in_7),Iso_k_95prc);
[th_g_in_80]=find_final_th_g(7*24,1,Iso_th(in_7),Iso_thp(in_7),Iso_k_80prc);
[th_g_in_65]=find_final_th_g(7*24,1,Iso_th(in_7),Iso_thp(in_7),Iso_k_65prc);
[th_g_in_5]=find_final_th_g(7*24,1,Iso_th(in_1),Iso_thp(in_1),Iso_k_5prc);
[th_g_in_20]=find_final_th_g(7*24,1,Iso_th(in_1),Iso_thp(in_1),Iso_k_20prc);
[th_g_in_35]=find_final_th_g(7*24,1,Iso_th(in_1),Iso_thp(in_1),Iso_k_35prc);

% mark out shaded area for intermidiate
inBetween_inte_90=[th_g_in_5, fliplr(th_g_in_95)];
inBetween_inte_60=[th_g_in_20, fliplr(th_g_in_80)];
inBetween_inte_30=[th_g_in_35, fliplr(th_g_in_65)];

% calculate upper and lower bound of peripheral
[th_g_pe_95]=find_final_th_g(7*24,1,Iso_th(pe_7),Iso_thp(pe_7),Iso_k_95prc);
[th_g_pe_80]=find_final_th_g(7*24,1,Iso_th(pe_7),Iso_thp(pe_7),Iso_k_80prc);
[th_g_pe_65]=find_final_th_g(7*24,1,Iso_th(pe_7),Iso_thp(pe_7),Iso_k_65prc);
[th_g_pe_5]=find_final_th_g(7*24,1,Iso_th(pe_1),Iso_thp(pe_1),Iso_k_5prc);
[th_g_pe_20]=find_final_th_g(7*24,1,Iso_th(pe_1),Iso_thp(pe_1),Iso_k_20prc);
[th_g_pe_35]=find_final_th_g(7*24,1,Iso_th(pe_1),Iso_thp(pe_1),Iso_k_35prc);

% mark out shaded area for peripheral
inBetween_peri_90=[th_g_pe_5, fliplr(th_g_pe_95)];
inBetween_peri_60=[th_g_pe_20, fliplr(th_g_pe_80)];
inBetween_peri_30=[th_g_pe_35, fliplr(th_g_pe_65)];


figure(2)
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis([0 7 1 1.65]);
hold on
fill(x_fill, inBetween_apex_90, 'r','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_apex_60, 'r','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_apex_30, 'r','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[Iso_thg(ap_1),Iso_thg(ap_3),Iso_thg(ap_7)],'r','filled');
fill(x_fill, inBetween_inte_90, 'b','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_inte_60, 'b','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_inte_30, 'b','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[Iso_thg(in_1),Iso_thg(in_3),Iso_thg(in_7)],'b','filled');
fill(x_fill, inBetween_peri_90, 'k','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_peri_60, 'k','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_peri_30, 'k','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[Iso_thg(pe_1),Iso_thg(pe_3),Iso_thg(pe_7)],'k','filled');

%% Plot of th_g over time with Abaqus simulation
%txt file each row is 95,80,65,(50 only for apex),35,20,5 percentile of th_g over time in Days
%
Iso_apex=load('Iso_apex_simulation.txt');
Iso_inte=load('Iso_inte_simulation.txt');
Iso_peri=load('Iso_peri_simulation.txt');

x_fill=[[0:7] fliplr([0:7])];
inBetween_apex_90=[Iso_apex(7,:), fliplr(Iso_apex(1,:))];
inBetween_apex_60=[Iso_apex(6,:), fliplr(Iso_apex(2,:))];
inBetween_apex_30=[Iso_apex(5,:), fliplr(Iso_apex(3,:))];

inBetween_inte_90=[Iso_inte(6,:), fliplr(Iso_inte(1,:))];
inBetween_inte_60=[Iso_inte(5,:), fliplr(Iso_inte(2,:))];
inBetween_inte_30=[Iso_inte(4,:), fliplr(Iso_inte(3,:))];

inBetween_peri_90=[Iso_peri(6,:), fliplr(Iso_peri(1,:))];
inBetween_peri_60=[Iso_peri(5,:), fliplr(Iso_peri(2,:))];
inBetween_peri_30=[Iso_peri(4,:), fliplr(Iso_peri(3,:))];

figure(3)
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis([0 7 1 1.65]);
hold on
fill(x_fill, inBetween_apex_90, 'r','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_apex_60, 'r','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_apex_30, 'r','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[Iso_thg(ap_1),Iso_thg(ap_3),Iso_thg(ap_7)],'r','filled');
fill(x_fill, inBetween_inte_90, 'b','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_inte_60, 'b','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_inte_30, 'b','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[Iso_thg(in_1),Iso_thg(in_3),Iso_thg(in_7)],'b','filled');
fill(x_fill, inBetween_peri_90, 'k','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_peri_60, 'k','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_peri_30, 'k','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[Iso_thg(pe_1),Iso_thg(pe_3),Iso_thg(pe_7)],'k','filled');