import os
os.environ["OMP_NUM_THREADS"] = "1"
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import numpy as np
import pandas as pd
import numpy.core.multiarray
import fastcluster
from multiprocessing import Pool,cpu_count
from functools import partial
from contextlib import closing
from joblib import Parallel, delayed

def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i
def flatx(x):
    return tuple(sorted(flatten(x)))

def dist(R,method):
	N = R.shape[0]
	d = R[np.triu_indices(N,1)]


	if method=='average':
		out = fastcluster.average(d)
	if method=='complete':
		out = fastcluster.complete(d)
	if method=='single':
		out = fastcluster.single(d)

	outI = out.astype(int)

	dend = {i:(i,) for i in xrange(N)}

	for i in xrange(len(outI)):
		dend[i+N] = (dend[outI[i][0]],dend[outI[i][1]])


	for i in xrange(N):
		dend.pop(i,None)

	dend = [(flatx(a),flatx(b)) for a,b in dend.values()]

	dend ={flatx((a,b)):(np.array(a),np.array(b)) for a,b in dend}
	return dend
   
def singPV(LV,Rb,method,gen):
    
    if method=='average':
        rxy = dict(map(lambda (k,(a,b)): (k,Rb[a][:,b].mean()),LV.items()))
    if method=='complete':
        rxy = dict(map(lambda (k,(a,b)): (k,Rb[a][:,b].max()),LV.items()))
    if method=='single':
        rxy = dict(map(lambda (k,(a,b)): (k,Rb[a][:,b].min()),LV.items()))

    PV = map(lambda (k,q):int(rxy[q]<rxy[k]), gen.items())
    
    return PV

def SingleBoot(X,LV,nan,method,gen,sel):
    
    Xb = X[:,sel]
    if nan==False:
        Rb = 1.-np.corrcoef(Xb)
    else:
        Rb = 1.-np.array(pd.DataFrame(Xb.T).corr())
        
    return singPV(LV,Rb,method,gen)

    
def BootDist(X,LV,method,Nt=1000,nan=False,ncpu=1):

	gen = {tuple(s):k for k in LV for s in LV[k] if len(s)>1}

	f = partial(SingleBoot,X,LV,nan,method,gen)

	sels = np.random.choice(range(X.shape[1]),replace=True,size=(Nt,X.shape[1]))


	PV = Parallel(n_jobs=ncpu, backend="threading")(delayed(f)(sel) for sel in sels)

	return sorted(zip(np.array(PV).mean(axis=0),gen.keys()))
    
def FDR(PV,alpha,N):
	p = np.array(zip(*PV)[0])
	thr = np.arange(1,len(PV)+1)*alpha/len(PV)

	sel = np.where(p<thr)[0]
	if len(sel)==0: 
		thr=-1.
	else:
		thr = thr[sel][-1]

	L = [range(N)]+[PV[i][1] for i in np.where(p<=thr)[0]]
	L = map(tuple,sorted(L,key=len,reverse=True))
	
	PV = {c:p for p,c in PV}
	PV[tuple(range(N))] = np.nan
	return L,PV
	
	
def Find_ValidatedCluster(X,Nt=1000,alpha=0.05,nan=False,ncpu=cpu_count(),method='average'):
	if nan==False:
		R = 1 - np.corrcoef(X)
	else:
		R = 1. -np.array(pd.DataFrame(X.T).corr()) 
			
	LV = dist(R,method)

	PV = BootDist(X,LV,method,Nt,nan,ncpu)
	
	L,PV = FDR(PV,alpha,X.shape[0])

	return L,PV,LV
