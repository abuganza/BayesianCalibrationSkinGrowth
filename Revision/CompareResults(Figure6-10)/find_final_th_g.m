function [th_g_vec]=find_final_th_g(t,dt,th,th_p,k)
t_vec=[0:dt:t];
th_g_dot_vec=zeros(1,length(t_vec));
th_g_vec=ones(1,length(t_vec));

for i=2:length(t_vec)
    th_g_dot_vec(i-1)=k*(th*th_p/th_g_vec(i-1)-th_p);
    th_g_vec(i)=th_g_vec(i-1)+dt*th_g_dot_vec(i-1);
end