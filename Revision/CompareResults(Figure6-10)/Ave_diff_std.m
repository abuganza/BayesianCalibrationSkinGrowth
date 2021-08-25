function [ave_diff,output] = Ave_diff_std(k,tdata,dt,th_data,thp_ave,thg_data)
%% this function calculate the ave.% of erors given a set of data
% input k is the growth parameter [hr^-1]
% input tdata,th_data,thg_data are vectors of input data
% input thp_ave is scaler: average theta_p
% output ave_diff is scaler:average of errors between th_g data and predict
% output output is vector of predicted th_g, same size as thg_data

num=length(tdata);
diff=zeros(1,num);  % percentage diff compare to data, same size as th_data
output=zeros(1,num);
for counter=1:num
    th_g_vec=find_final_th_g(tdata(counter),dt,th_data(counter),thp_ave,k);
    diff(counter)=abs(thg_data(counter)-th_g_vec(end))/thg_data(counter);
    output(1,counter)=th_g_vec(end);
end
ave_diff=mean(diff);
end

