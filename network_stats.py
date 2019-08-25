import micro_nw as nw
import numpy as np
from sklearn.linear_model import Lasso, LassoCV
import pandas as pd
import sys
from scipy import stats

folde = sys.argv[1]
fo1 = open(folde+'/disease_list','w')
fo1.write(sys.argv[2])
dis =[]
dis.append(sys.argv[2])
if len(sys.argv)>3:
	dis.append(sys.argv[3])
	fo1.write('\t'+sys.argv[3])
if len(sys.argv)>4:
	dis.append(sys.argv[4])
	fo1.write('\t'+sys.argv[4])
fo1.close()


alp = 0.0001
reps = 1000



D = pd.read_csv(folde+'/healthy.csv').as_matrix().astype('float64')
max_degree ,components, degrees, GG = nw.max_deg_stat(D,alp,reps)






for d in dis:
	D= pd.read_csv(folde+'/'+d).as_matrix().astype('float64')
	md ,c, d, g = nw.max_deg_stat(D,alp,reps)
	max_degree = np.concatenate((max_degree,md),axis=0)
	components = np.concatenate((components,c),axis=0)
	degrees = np.concatenate((degrees,d),axis=0)
	GG=GG+g

np.savetxt(folde+'/degrees.csv',degrees,delimiter=',')

R,P = nw.rank_corr(degrees)
np.savetxt(folde+'/degree_r.csv',R,delimiter=',')
np.savetxt(folde+'/degree_p.csv',P,delimiter=',')
S=nw.graph_similarity(GG)
np.savetxt(folde+'/graph_similar.csv',S,delimiter=',')