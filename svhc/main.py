import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import pandas as pd
import numpy.core.multiarray
import fastcluster
from multiprocessing import Pool,cpu_count
from functools import partial
from contextlib import closing

def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i
def flatx(x):
    return tuple(sorted(flatten(x)))

def AVdist(R):
    N = R.shape[0]
    d = R[np.triu_indices(N,1)]
    out = fastcluster.average(d)

    outI = out.astype(int)

    dend = {i:(i,) for i in xrange(N)}

    for i in xrange(len(outI)):
        dend[i+N] = (dend[outI[i][0]],dend[outI[i][1]])


    for i in xrange(N):
        dend.pop(i,None)

    dend = [(flatx(a),flatx(b)) for a,b in dend.values()]

    dend ={flatx((a,b)):(np.array(a),np.array(b)) for a,b in dend}
    return dend
   
def SingleBoot(X,LV,nan,seed=None):
	local_state = np.random.RandomState(seed)
	sel = local_state.choice(range(X.shape[1]),replace=True,size=X.shape[1])
	Xb = X[:,sel]
	if nan==False:
		Rb = 1.-np.corrcoef(Xb)
	else:
		Rb = 1.-np.array(pd.DataFrame(Xb.T).corr())
	return singPV(LV,Rb)

    
def BootDist(X,LV,Nt=1000,nan=False,ncpu=1):

	f = partial(SingleBoot,X,LV,nan)
	seeds = np.random.randint(2**32,size=Nt)
	if ncpu==1:
		PV = map(f,seeds)
	else:
		with closing(Pool(processes=ncpu)) as p:
			PV = list(p.imap_unordered(f,seeds))

	return PV

def singPV(LV,Rb):

    vl = [(tuple(a),tuple(b)) for a,b in LV.values()]
    rxy = dict(zip(vl,map(lambda (a,b): Rb[a][:,b].mean(),LV.values())))

    PV = []
    for (a,b) in LV.values():
        a,b = tuple(a),tuple(b)
        t = b
        if LV.has_key(t):
            c,d = map(tuple,LV[t])
            PV.append((((rxy[(c,d)]-rxy[(a,b)])>0).sum(),t))
        t = a
        if LV.has_key(t):
            c,d = map(tuple,LV[t])
            PV.append((((rxy[(c,d)]-rxy[(a,b)])>0).sum(),t))
    return zip(*PV)[0]

def Validate(PV,LV,Nt,N,alpha):
    PV = (np.array(PV).sum(axis=0))/float(Nt)

    nod = [tuple(e) for (a,b) in LV.values() for e in [b,a] if LV.has_key(tuple(e))]

    PV = zip(PV,nod)

    PV = sorted(PV)

    L,PV = FDR(PV,alpha,N)
    
    return L,PV


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

def Find_ValidatedCluster(X,Nt=1000,alpha=0.05,nan=False,ncpu=cpu_count()):
	if nan==False:
		R = 1 - np.corrcoef(X)
	else:
		R = 1. -np.array(pd.DataFrame(X.T).corr()) 

	LV = AVdist(R)

	PV = BootDist(X,LV,Nt,nan,ncpu)
	
	L,PV = Validate(PV,LV,Nt,X.shape[0],alpha)

	
	return L,PV,LV
