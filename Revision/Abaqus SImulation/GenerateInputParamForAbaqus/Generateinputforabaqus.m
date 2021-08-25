clc
clear all;
close all;

%% For Iso Cases
possion=0.349;
mu_Iso_array=0.74*[0.9 1 1.1];
lam_Iso_array=0.544*[0.9 1 1.1];
pymc3_Iso=load('pymc3_Iso_output_v2.txt');
k_Iso=exp(pymc3_Iso)/25;
k_Iso_array=24*prctile(k_Iso,[5 15 25 35 45 55 65 75 85 95]); % Translate from hr^-1 to day^-1
%% Generate matrix of input
Iso_input=zeros([30,3]);
Iso_input(:,3)=repmat(k_Iso_array',[3,1]);
lam_repeated=repelem(lam_Iso_array,10);
mu_repeated=repelem(mu_Iso_array,10);
Iso_input(:,1)=lam_repeated';
Iso_input(:,2)=mu_repeated';

%% For Aniso Cases
pymc3_G1=load('pymc3_G1_output_v2.txt');
pymc3_G2=load('pymc3_G2_output_v2.txt');
k_G1=exp(pymc3_G1)/50;
k_G2=exp(pymc3_G2)/0.5;
mu_array=0.05*[0.9 1 1.1];
k_array=4.8*mu_array;
k_G1_array=24*prctile(k_G1,[10 30 50 70 90]);
k_G2_array=24*prctile(k_G2,[10 30 50 70 90]);
%% Generate Matrix of input
Anis_input=zeros([75,4]);
mu_repeated=repelem(mu_array,25);
k_repeated=repelem(k_array,25);
k_G1_repeated=repmat(repelem(k_G1_array,5),[1,3]);
k_G2_repeated=repmat(k_G2_array,[1,15]);
Anis_input(:,1)=mu_repeated';
Anis_input(:,2)=k_repeated';
Anis_input(:,3)=k_G1_repeated';
Anis_input(:,4)=k_G2_repeated';

%%
writematrix(Iso_input,'Iso_lam_mu_k','Delimiter','space')
writematrix(Anis_input,'Anis_mu_k_G1_G2','Delimiter','space')