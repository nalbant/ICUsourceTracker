import numpy as np
from sklearn.linear_model import Lasso, LassoCV
import pandas as pd
import sys
from scipy import stats
import scipy.sparse.csgraph as gr
from scipy import sparse

def max_degree(G,le):
	deg=0
	for i in range(le):
		temp_deg = (G[i,:]!=0).sum() + (G[:,i]!=0).sum()
		if temp_deg > deg:
			deg = temp_deg
	return deg

def degree_list(G,le):
	degrees=[]
	for i in range(le):
		degrees.append((G[i,:]!=0).sum() + (G[:,i]!=0).sum())
	return np.array(degrees)


def adj_matrix(x,alp):
	n_otu = x.shape[1]
	G=np.zeros((n_otu,n_otu))
	for i in range(n_otu):
		A=x.copy()
		y = A[:,i].copy()
		A[:,i]=0.0
		model = Lasso(alpha=alp, fit_intercept=False, tol=0.0000001).fit(A, y)
		G[i] = coef = model.coef_
	return G



def max_deg_stat(x,alp,reps):
	n_otu = x.shape[1]
	n_samp = x.shape[0]
	boot = 20#int(n_samp/2)
	mdeg = np.zeros(reps)
	components = np.zeros(reps)
	degrees= np.zeros((reps,n_otu))
	GG=[]

	for i in range(reps):
		idx = np.random.permutation(n_samp)[:boot]
		G = adj_matrix(x[idx,:],alp)
		GG.append(G)
		mdeg[i] = max_degree(G,n_otu)
		components[i] = gr.connected_components(G, return_labels=False)
		degrees[i] = degree_list(G,n_otu)
	return mdeg ,components, degrees, GG


def network_ensemble(x,alp,reps):
	n_otu = x.shape[1]
	n_samp = x.shape[0]
	boot = 20#int(n_samp/2)
	GG=[]

	for i in range(reps):
		idx = np.random.permutation(n_samp)[:boot]
		G = sparse.csr_matrix(adj_matrix(x[idx,:],alp))
		GG.append(G)
	return GG

def network_incremental(x,alp,reps):
	n_samp = x.shape[0]
	if n_samp<100:
		boot= list(range(2,n_samp))
	else:
		boot= np.linspace(2,n_samp-1,num=100,dtype=int)
	GG=[]
	for bb in boot:
		for i in range(reps):
			idx = np.random.permutation(n_samp)[:bb]
			G = sparse.csr_matrix(adj_matrix(x[idx,:],alp))
			GG.append(G)
	G = sparse.csr_matrix(adj_matrix(x,alp))
	GG.append(G)
	return GG

def rank_corr(A):
	n = A.shape[0]
	R = np.zeros((n,n))
	P = np.zeros((n,n))
	for i in range(n-1):
		R[i,i]=0.5
		for j in range(i+1,n):
			r,p = stats.spearmanr(A[i],A[j])
			R[i,j] = r
			P[i,j] = p
	R=R+R.T
	R[-1,-1]=1
	P=P+P.T
	return R,P

def graph_similarity(GG):
	n = len(GG)
	S = np.zeros((n,n))
	for i in range(n-1):
		S[i,i]=0.5
		for j in range(i+1,n):
			S[i,j] = (GG[i]*GG[j]!=0).sum()/((GG[i]+GG[j])!=0).sum()
	S=S+S.T
	S[-1,-1]=1
	return S

# def zipfs_exponent(degrees):
# 	n = len(degrees)
# 	zfe = np.zeros(n)
# 	for i in range(n):
# 		ss= np.sort(degrees[i])[::-1]
# 		ss = np.log(ss[ss>0])
