import numpy as np
import igraph as ig
import sys
import pandas as pd
import hclustval
import numpy.core.multiarray
import fastcluster


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
    
def BootDist(X,Nt=1000,nan=False):

    Rb = []
    for i in xrange(Nt):
        sel = np.random.choice(range(X.shape[1]),replace=True,size=X.shape[1])
        Xb = X[:,sel]
        if nan==False:
            Rb.append(1.-np.corrcoef(Xb))
        else:
            Rb.append(1.-np.array(pd.DataFrame(Xb.T).corr()))

    return Rb
    

def HclustValDist(X,Nt=1000,alpha=0.05,nan=False):
	if nan==False:
		R = 1 - np.corrcoef(X)
	else:
		R = 1. -np.array(pd.DataFrame(X.T).corr()) 

	LV = AVdist(R)

	Rb = BootDist(X,Nt,nan)

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

	p = np.array(zip(*PV)[0])

	thr = np.arange(1,len(PV)+1)*alpha/len(PV)


	sel = np.where(p<thr)[0]
	if len(sel)==0: 
		thr=-1.
	else:
		thr = thr[sel][-1]

	N = X.shape[0]
	L = [range(N)]+[PV[i][1] for i in np.where(p<=thr)[0]]
	L = map(tuple,sorted(L,key=len,reverse=True))
	
	return L,LV
