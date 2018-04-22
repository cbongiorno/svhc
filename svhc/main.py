import numpy as np
import pandas as pd
import numpy.core.multiarray
import fastcluster
from multiprocessing import Pool
from functools import partial


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
   
def SingleBoot(X,nan,seed=None):
	local_state = np.random.RandomState(seed)
	sel = local_state.choice(range(X.shape[1]),replace=True,size=X.shape[1])
	Xb = X[:,sel]
	if nan==False:
		return 1.-np.corrcoef(Xb)
	else:
		return 1.-np.array(pd.DataFrame(Xb.T).corr())

    
def BootDist(X,Nt=1000,nan=False,ncpu=1):

	f = partial(SingleBoot,X,nan)
	seeds = np.random.randint(2**32,size=Nt)
	if ncpu==1:
		Rb = map(f,seeds)
	else:
		p = Pool(processes=ncpu)
		Rb = p.map(f,seeds)
		p.close()

	return Rb
    
def get_Pvalues(Rb,LV,N,alpha):
	vl = [(tuple(a),tuple(b)) for a,b in LV.values()]
	rxy = dict(zip(vl,map(lambda (a,b): np.array(map(lambda x:x[a][:,b].mean(),Rb)),LV.values())))

	Nt = float(len(Rb))
	PV = []
	for (a,b) in LV.values():
		a,b = tuple(a),tuple(b)
		t = b
		if LV.has_key(t):
			c,d = map(tuple,LV[t])
			PV.append((((rxy[(c,d)]-rxy[(a,b)])>0).sum()/Nt,t))
		t = a
		if LV.has_key(t):
			c,d = map(tuple,LV[t])
			PV.append((((rxy[(c,d)]-rxy[(a,b)])>0).sum()/Nt,t))
		
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

def Find_ValidatedCluster(X,Nt=1000,alpha=0.05,nan=False,ncpu=1):
	if nan==False:
		R = 1 - np.corrcoef(X)
	else:
		R = 1. -np.array(pd.DataFrame(X.T).corr()) 

	LV = AVdist(R)

	Rb = BootDist(X,Nt,nan,ncpu)
	
	L,PV = get_Pvalues(Rb,LV,X.shape[0],alpha)


	
	return L,PV,LV
