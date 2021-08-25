clc
clear all;
close all;

%% Load Noise data
Iso_noise=load('Iso_post_noise.txt');
G1_noise=load('G1_post_noise.txt');
G2_noise=load('G2_post_noise.txt');

%% Histogram plot
figure(1)
hold on
histogram(Iso_noise,40,'Normalization','probability','FaceColor',[0.4 0.4 0.4]);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
xlabel('Posterior Noise Isotropic')
ylabel('Probability')
hold off

figure(2)
hold on
histogram(G1_noise,40,'Normalization','probability','FaceColor',[0.4 0.4 0.4]);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
xlabel('Posterior Noise G_1')
ylabel('Probability')
hold off

figure(3)
hold on
histogram(G2_noise,40,'Normalization','probability','FaceColor',[0.4 0.4 0.4]);
box on
grid minor
ax=gca;
ax.FontSize = 18;
ax.FontName='Arial';
xlabel('Posterior Noise G_2')
ylabel('Probability')
hold off