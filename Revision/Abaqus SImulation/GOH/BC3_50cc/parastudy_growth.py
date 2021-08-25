import numpy as np
import os
import shutil

# load parameter list generated by LHS
# mu, n, k, m

paralist = np.loadtxt('Anis_mu_k_G1_G2.txt')
nu = 0.4
tmp=451
for ParamVal in paralist:
	parfile = open('param.param','w')
	parfile.write('*PARAMETER\n')
	parfile.write('p=0.03\n')
	parfile.write('t_init=0.001\n')
	parfile.write('t_max=0.001\n')
	parfile.write('k= ')
	parfile.write('%.3f\n'% ParamVal[1])
	parfile.write('mu= ')
	parfile.write('%.3f\n'% ParamVal[0])
	parfile.write('kappa=0.0498\n')
	parfile.write('k1= 4.97\n')
	parfile.write('k2=2.88\n')
	parfile.write('f0_1=1.\n')
	parfile.write('f0_2=0.\n')
	parfile.write('f0_3=0.\n')
	parfile.write('xn0_1=0.\n')
	parfile.write('xn0_2=0.\n')
	parfile.write('xn0_3=1.\n')
	parfile.write('kkg1= ')
	parfile.write('%.3f\n'% ParamVal[2])
	parfile.write('kkg2= ')
	parfile.write('%.3f\n'% ParamVal[3])
	parfile.write('mg1=0.\n')	
	parfile.write('mg2=0.\n')
	parfile.write('ng1=0.\n')	
	parfile.write('ng2=0.\n')
	parfile.write('tcr1=1.1254\n')		
	parfile.write('tcr2=1.1175\n')	
	parfile.close()
	os.system("abaqus job=JobName input=GOH_BC3 user=GOH_50cc ask_delete=OFF cpus=24 interactive")
	os.system("abaqus viewer noGUI=readODB.py")

	os.rename("JobName.sta", "job" + str(tmp)+ ".sta")
	os.rename("JobName.odb", "job" + str(tmp)+ ".odb")
	os.rename("JobName.pes", "job" + str(tmp)+ ".pes")
	os.rename("JobName1.txt", "job" + str(tmp)+ "lam1g.txt")
	os.rename("JobName2.txt", "job" + str(tmp)+ "lam1.txt")
	os.rename("JobName3.txt", "job" + str(tmp)+ "lam2g.txt")
	os.rename("JobName4.txt", "job" + str(tmp)+ "lam2.txt")
	os.rename("JobName5.txt", "job" + str(tmp)+ "xcoord.txt")
	os.rename("JobName6.txt", "job" + str(tmp)+ "ycoord.txt")
	os.rename("JobName7.txt", "job" + str(tmp)+ "zcoord.txt")
	tmp += 1
	#os.rename("JobName.odb", "job_"+ str(int(ParamVal[0]*1000.))+"_"+ str(int(ParamVal[1]*1000.))+"_"+ str(int(ParamVal[2]*1000.)) +"_"+ str(int(ParamVal[3]*1000.)) + ".odb")
