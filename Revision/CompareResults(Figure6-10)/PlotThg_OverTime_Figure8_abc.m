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
histogram(k_G2,40,'Normalization','probability','FaceColor',[0.4 0.4 0.4]);
G2_k_5prc=prctile(k_G2,5);
G2_k_20prc=prctile(k_G2,20);
G2_k_35prc=prctile(k_G2,35);
G2_k_50prc=prctile(k_G2,50);
G2_k_65prc=prctile(k_G2,65);
G2_k_80prc=prctile(k_G2,80);
G2_k_95prc=prctile(k_G2,95);

%% plot of G2 th_g over time
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
ap_1=30;
ap_3=59;
ap_7=106; %246 200 196

in_1=19;
in_3=51; 
in_7=102;

pe_1=9;
pe_3=55;
pe_7=219;

% calculate upper and lower bound of apex
[th_g_ap_95]=find_final_th_g(7*24,1,G2_th(ap_7),G2_thp(ap_7),G2_k_95prc);
[th_g_ap_80]=find_final_th_g(7*24,1,G2_th(ap_7),G2_thp(ap_7),G2_k_80prc);
[th_g_ap_65]=find_final_th_g(7*24,1,G2_th(ap_7),G2_thp(ap_7),G2_k_65prc);
[th_g_ap_5]=find_final_th_g(7*24,1,G2_th(ap_3),G2_thp(ap_3),G2_k_5prc);
[th_g_ap_20]=find_final_th_g(7*24,1,G2_th(ap_3),G2_thp(ap_3),G2_k_20prc);
[th_g_ap_35]=find_final_th_g(7*24,1,G2_th(ap_3),G2_thp(ap_3),G2_k_35prc);

% mark out shaded area for apex
x_fill=[[0:7*24]/24 fliplr([0:7*24]/24)];
inBetween_apex_90=[th_g_ap_5, fliplr(th_g_ap_95)];
inBetween_apex_60=[th_g_ap_20, fliplr(th_g_ap_80)];
inBetween_apex_30=[th_g_ap_35, fliplr(th_g_ap_65)];

% calculate upper and lower bound of intermidiate
[th_g_in_95]=find_final_th_g(7*24,1,G2_th(in_7),G2_thp(in_7),G2_k_95prc);
[th_g_in_80]=find_final_th_g(7*24,1,G2_th(in_7),G2_thp(in_7),G2_k_80prc);
[th_g_in_65]=find_final_th_g(7*24,1,G2_th(in_7),G2_thp(in_7),G2_k_65prc);
[th_g_in_5]=find_final_th_g(7*24,1,G2_th(in_1),G2_thp(in_1),G2_k_5prc);
[th_g_in_20]=find_final_th_g(7*24,1,G2_th(in_1),G2_thp(in_1),G2_k_20prc);
[th_g_in_35]=find_final_th_g(7*24,1,G2_th(in_1),G2_thp(in_1),G2_k_35prc);

% mark out shaded area for intermidiate
inBetween_inte_90=[th_g_in_5, fliplr(th_g_in_95)];
inBetween_inte_60=[th_g_in_20, fliplr(th_g_in_80)];
inBetween_inte_30=[th_g_in_35, fliplr(th_g_in_65)];

% calculate upper and lower bound of peripheral
[th_g_pe_95]=find_final_th_g(7*24,1,G2_th(pe_7),G2_thp(pe_7),G2_k_95prc);
[th_g_pe_80]=find_final_th_g(7*24,1,G2_th(pe_7),G2_thp(pe_7),G2_k_80prc);
[th_g_pe_65]=find_final_th_g(7*24,1,G2_th(pe_7),G2_thp(pe_7),G2_k_65prc);
[th_g_pe_5]=find_final_th_g(7*24,1,G2_th(pe_1),G2_thp(pe_1),G2_k_5prc);
[th_g_pe_20]=find_final_th_g(7*24,1,G2_th(pe_1),G2_thp(pe_1),G2_k_20prc);
[th_g_pe_35]=find_final_th_g(7*24,1,G2_th(pe_1),G2_thp(pe_1),G2_k_35prc);

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
axis([0 7 1 1.35]);
hold on
fill(x_fill, inBetween_apex_90, 'r','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_apex_60, 'r','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_apex_30, 'r','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[G2_thg(ap_1),G2_thg(ap_3),G2_thg(ap_7)],'r','filled');
fill(x_fill, inBetween_inte_90, 'b','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_inte_60, 'b','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_inte_30, 'b','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[G2_thg(in_1),G2_thg(in_3),G2_thg(in_7)],'b','filled');
fill(x_fill, inBetween_peri_90, 'k','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_peri_60, 'k','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_peri_30, 'k','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[G2_thg(pe_1),G2_thg(pe_3),G2_thg(pe_7)],'k','filled');
%% Plot of th_g over time with Abaqus simulation
%txt file each row is 95,80,65,(50 only for apex),35,20,5 percentile of th_g over time in Days
G2_apex=load('G2_apex_simulation.txt');
G2_inte=load('G2_inte_simulation.txt');
G2_peri=load('G2_peri_simulation.txt');

x_fill=[[0:7] fliplr([0:7])];
inBetween_apex_90=[G2_apex(7,:), fliplr(G2_apex(1,:))];
inBetween_apex_60=[G2_apex(6,:), fliplr(G2_apex(2,:))];
inBetween_apex_30=[G2_apex(5,:), fliplr(G2_apex(3,:))];

inBetween_inte_90=[G2_inte(6,:), fliplr(G2_inte(1,:))];
inBetween_inte_60=[G2_inte(5,:), fliplr(G2_inte(2,:))];
inBetween_inte_30=[G2_inte(4,:), fliplr(G2_inte(3,:))];

inBetween_peri_90=[G2_peri(6,:), fliplr(G2_peri(1,:))];
inBetween_peri_60=[G2_peri(5,:), fliplr(G2_peri(2,:))];
inBetween_peri_30=[G2_peri(4,:), fliplr(G2_peri(3,:))];

figure(3)
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
axis([0 7 1 1.35]);
hold on
fill(x_fill, inBetween_apex_90, 'r','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_apex_60, 'r','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_apex_30, 'r','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[G2_thg(ap_1),G2_thg(ap_3),G2_thg(ap_7)],'r','filled');
fill(x_fill, inBetween_inte_90, 'b','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_inte_60, 'b','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_inte_30, 'b','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[G2_thg(in_1),G2_thg(in_3),G2_thg(in_7)],'b','filled');
fill(x_fill, inBetween_peri_90, 'k','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_peri_60, 'k','FaceAlpha',.1,'EdgeColor','none');
fill(x_fill, inBetween_peri_30, 'k','FaceAlpha',.1,'EdgeColor','none');
scatter([1,3,7],[G2_thg(pe_1),G2_thg(pe_3),G2_thg(pe_7)],'k','filled');