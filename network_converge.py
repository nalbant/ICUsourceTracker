import micro_nw as nw
import numpy as np
import pandas as pd
import sys
import pickle


params = sys.argv[:-1]
folde = params[1]
idx = sys.argv[-1]
fo1 = open(folde+'/disease_list','w')
fo1.write(params[2])
dis =[]
dis.append(params[2])
if len(params)>3:
	dis.append(params[3])
	fo1.write('\t'+params[3])
if len(params)>4:
	dis.append(params[4])
	fo1.write('\t'+params[4])
fo1.close()


alp = 0.0001
reps = 5



D = pd.read_csv(folde+'/healthy.csv').as_matrix().astype('float64')
GG=[]
GG.append(nw.network_incremental(D,alp,reps))

for d in dis:
	D= pd.read_csv(folde+'/'+d).as_matrix().astype('float64')
	g = nw.network_incremental(D,alp,reps)
	GG.append(g)

fo = open(folde+'/'+idx+'_converge.pck','wb')
pickle.dump(GG,fo)
fo.close()

#fi = open(folde+'/'+idx+'.pck','rb')
#al = pickle.load(fi)
#fi.close()