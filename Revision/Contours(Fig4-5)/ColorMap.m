clc
clear all
close all
%% Load data
th_1h=load('th_1h.txt');
th_1d=load('th_1d.txt');
th_3d=load('th_3d.txt');
th_7d=load('th_7d.txt');

th_g_1h=load('th_g_1h.txt');
th_g_1d=load('th_g_1d.txt');
th_g_3d=load('th_g_3d.txt');
th_g_7d=load('th_g_7d.txt');

th_p_1h=load('th_p_1h.txt');
th_p_1d=load('th_p_1d.txt');
th_p_3d=load('th_p_3d.txt');
th_p_7d=load('th_p_7d.txt');

th_e_1h=load('th_e_1h.txt');
th_e_1d=load('th_e_1d.txt');
th_e_3d=load('th_e_3d.txt');
th_e_7d=load('th_e_7d.txt');

th_1h_G1=load('th_1h_G1.txt');
th_1d_G1=load('th_1d_G1.txt');
th_3d_G1=load('th_3d_G1.txt');
th_7d_G1=load('th_7d_G1.txt');

th_g_1h_G1=load('th_g_1h_G1.txt');
th_g_1d_G1=load('th_g_1d_G1.txt');
th_g_3d_G1=load('th_g_3d_G1.txt');
th_g_7d_G1=load('th_g_7d_G1.txt');

th_p_1h_G1=load('th_p_1h_G1.txt');
th_p_1d_G1=load('th_p_1d_G1.txt');
th_p_3d_G1=load('th_p_3d_G1.txt');
th_p_7d_G1=load('th_p_7d_G1.txt');

th_e_1h_G1=load('th_e_1h_G1.txt');
th_e_1d_G1=load('th_e_1d_G1.txt');
th_e_3d_G1=load('th_e_3d_G1.txt');
th_e_7d_G1=load('th_e_7d_G1.txt');

th_1h_G2=load('th_1h_G2.txt');
th_1d_G2=load('th_1d_G2.txt');
th_3d_G2=load('th_3d_G2.txt');
th_7d_G2=load('th_7d_G2.txt');

th_g_1h_G2=load('th_g_1h_G2.txt');
th_g_1d_G2=load('th_g_1d_G2.txt');
th_g_3d_G2=load('th_g_3d_G2.txt');
th_g_7d_G2=load('th_g_7d_G2.txt');

th_p_1h_G2=load('th_p_1h_G2.txt');
th_p_1d_G2=load('th_p_1d_G2.txt');
th_p_3d_G2=load('th_p_3d_G2.txt');
th_p_7d_G2=load('th_p_7d_G2.txt');

th_e_1h_G2=load('th_e_1h_G2.txt');
th_e_1d_G2=load('th_e_1d_G2.txt');
th_e_3d_G2=load('th_e_3d_G2.txt');
th_e_7d_G2=load('th_e_7d_G2.txt');
%% Figure 4
figure(1)
set(gca,'FontSize',30)
subplot(2,5,1);
image((th_1h+th_1d+th_3d+th_7d)/4,'CDataMapping','scaled');
%colorbar

subplot(2,5,2)
image(th_g_1h,'CDataMapping','scaled');
caxis([0.6 2]);

subplot(2,5,3);
image(th_g_1d,'CDataMapping','scaled');
caxis([0.6 2]);

subplot(2,5,4);
image(th_g_3d,'CDataMapping','scaled');
caxis([0.6 2]);

subplot(2,5,5);
image(th_g_7d,'CDataMapping','scaled');
caxis([0.6 2]);

subplot(2,5,6);
image((th_p_1h+th_p_1d+th_p_3d+th_p_7d)/4,'CDataMapping','scaled');
%colorbar

subplot(2,5,7)
image(th_e_1h,'CDataMapping','scaled');
caxis([0.6 2]);

subplot(2,5,8);
image(th_e_1d,'CDataMapping','scaled');
caxis([0.6 2]);

subplot(2,5,9);
image(th_e_3d,'CDataMapping','scaled');
caxis([0.6 2]);

subplot(2,5,10);
image(th_e_7d,'CDataMapping','scaled');
caxis([0.6 2]);
colormap jet
%colorbar

%% Figure (5)
figure(2)
subplot(4,5,1);
image((th_1h_G1+th_1d_G1+th_3d_G1+th_7d_G1)/4,'CDataMapping','scaled');
caxis([0.8 1.3]);
%colorbar
subplot(4,5,2);
image(th_g_1h_G1,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,3);
image(th_g_1d_G1,'CDataMapping','scaled');
caxis([0.8 1.6]);
%colorbar
subplot(4,5,4);
image(th_g_3d_G1,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,5);
image(th_g_7d_G1,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,6);
image((th_p_1h_G1+th_p_1d_G1+th_p_3d_G1+th_p_7d_G1)/4,'CDataMapping','scaled');
caxis([0.8 1.3]);
%colorbar
subplot(4,5,7);
image(th_g_1h_G2,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,8);
image(th_g_1d_G2,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,9);
image(th_g_3d_G2,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,10);
image(th_g_7d_G2,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,11);
image((th_1h_G2+th_1d_G2+th_3d_G2+th_7d_G2)/4,'CDataMapping','scaled');
caxis([0.8 1.3]);
%colorbar
subplot(4,5,12);
image(th_e_1h_G1,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,13);
image(th_e_1d_G1,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,14);
image(th_e_3d_G1,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,15);
image(th_e_7d_G1,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,16);
image((th_p_1h_G2+th_p_1d_G2+th_p_3d_G2+th_p_7d_G2)/4,'CDataMapping','scaled');
caxis([0.8 1.3]);
%colorbar
subplot(4,5,17);
image(th_e_1h_G2,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,18);
image(th_e_1d_G2,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,19);
image(th_e_3d_G2,'CDataMapping','scaled');
caxis([0.8 1.6]);

subplot(4,5,20);
image(th_e_7d_G2,'CDataMapping','scaled');
caxis([0.8 1.6]);
colormap jet