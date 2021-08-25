
# coding: utf-8

# In[1]:


import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd 

import pymc3 as pm
import theano as T
from theano import function, shared, tensor as tt


import pandas as pd
import arviz as az

import matplotlib.pyplot as plt 
import seaborn as sns 
sns.set()
from graphviz import Digraph
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


data_day = pd.read_csv('Iso_Data_Filtered.csv')
data_day.dropna()
th = data_day.theta.values
thg = data_day.theta_g.values
time =data_day.time.values
#th.shape[0]


# In[3]:


#thg.shape[0]


# In[8]:


# for test run, only use first 5 data points
#th=th[0:5]
#thg=thg[0:5]
#time=time[0:5]

# In[11]:


"""
This function solve for the predicted th_g given a certain value of k.
output is an array same size at the thg
"""
def Solver(k):
  num_pts=np.size(th)
  k = tt.exp(k) / 25
  n=np.size(th);
  dt=1;
  predict=tt.zeros(num_pts)
  for counter in range(0, num_pts):
    temp_th=th[counter]
    temp_th_p=1.198
    temp_time=time[counter]
    th_g_pred=1;
    th_g_pred_s=shared(th_g_pred)

    for t_step in np.arange(0.,temp_time,dt):
      th_g_dot=k*(temp_th*temp_th_p/th_g_pred_s-temp_th_p)
      th_g_pred_s=th_g_pred_s+th_g_dot*dt;
    predict=tt.set_subtensor(predict[counter:counter+1], th_g_pred_s)

  return predict


# In[ ]:



model_C = pm.Model()
alpha1 =3.
beta1 = 0.05
alpha2 = 1.0
# define the distribution 
with model_C:
    sigma2s = pm.InverseGamma('sigma2s',alpha=alpha1,beta=beta1,shape=1)
    sigma2 = pm.Deterministic('sigma2',tt.tile(sigma2s,th.shape[0]))
    gamma2 = pm.Exponential(name='gamma2',lam=alpha2)
    ln_k_guess = pm.Normal(name='ln_k_guess',mu =0,sigma=tt.sqrt(gamma2),shape=1)
    y_mean = pm.Deterministic('y_mean', Solver(ln_k_guess))
    y = pm.Normal(name='y',mu=y_mean,sigma=tt.sqrt(sigma2),observed=thg)


# In[12]:


with model_C:
    mcmc_res_C = pm.sample(draws = 5000, step=pm.NUTS())


#_=pm.plot_posterior(mcmc_res_C, var_names=['ln_k_guess'])


# In[ ]:


#_=pm.plot_autocorr(mcmc_res_C, var_names=['ln_k_guess'], combined=True)


# In[ ]:


thinned_mcmc_res_C = mcmc_res_1[1000::10]
ppsamples_C = pm.sample_posterior_predictive(samples=1600,\
                                           model=model_C, 
                                           trace=mcmc_res_C_thin)['y']

#post_sample=pm.sampling.sample(model=model_C)

print(psamples_C[:]['ln_k_guess'])
np.savetxt('IsoKBayesian.txt',post_sample[:]['ln_k_guess'])
