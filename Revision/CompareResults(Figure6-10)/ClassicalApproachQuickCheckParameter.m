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

%% Use Classical Approach, calculate minimum ave error with different k and plot

%% Iso
num=40;
dt=0.1;
diff_ave=zeros(1,num);
k=0.001*[1:num]+0.02;
for ii =1:num
    [diff_ave(ii),temp]=Ave_diff_std(k(ii),Iso_time,dt,Iso_th,Iso_thp(1),Iso_thg);
end

figure(1)
plot(k,diff_ave);
xlabel('k')
ylabel('Iso Prediction Average Difference [%]')
grid minor
ax=gca;
ax.FontSize = 12;
%% G1
num=100;
dt=0.1;
diff_ave=zeros(1,num);
k=0.001*[1:num]+0.01;
for ii =1:num
    [diff_ave(ii),temp]=Ave_diff_std(k(ii),G1_time,dt,G1_th,G1_thp(1),G1_thg);
end

figure(2)
plot(k,diff_ave);
xlabel('k')
ylabel('G1 Prediction Average Difference [%]')
grid minor
ax=gca;
ax.FontSize = 12;

%% G2
num=60;
dt=0.1;
diff_ave=zeros(1,num);
k=0.001*[1:num]+0.01;
for ii =1:num
    [diff_ave(ii),temp]=Ave_diff_std(k(ii),G2_time,dt,G2_th,G2_thp(1),G2_thg);
end

figure(3)
plot(k,diff_ave);
xlabel('k')
ylabel('G2 Prediction Average Difference [%]')
grid minor
ax=gca;
ax.FontSize = 12;