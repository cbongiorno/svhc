import numpy as np
import igraph as ig

def ToMemb(nod,N):
    h = np.array([-1]*N)

    nod = sorted(nod,key=lambda x:len(x),reverse=True)
    Hn = []
    sel = {}
    while True:
        c = 0
        for i in range(len(nod)):
            if i in sel: continue
            if not all(h[np.array(nod[i])]==-1): continue

            h[np.array(nod[i])]=c
            sel[i] = None
            c+=1

        h[h==-1] = np.arange(h.max()+1,(h==-1).sum()+h.max()+1)
        Hn.append(h)
        h = np.array([-1]*N)
        if len(sel)==len(nod):
            break

    Hn.append(np.arange(N))    
    return Hn

def getPvalue(a,b):
    df =  ((a-b)>0).sum()
    if df>len(a)/2.:
        return (len(a)-df)/float(len(a))
    else:
        return df/float(len(a))

def srt(x):
	return tuple(sorted(x))

def AverageLink(R):
	N = R.shape[0]
	LV = {}

	nod =  [np.array([i]) for i in xrange(N)]

	for lv in xrange(N-1):
		Y = []
		for i in xrange(len(nod)):
			for j in xrange(i):
				Y.append((R[nod[i]][:,nod[j]].mean(),i,j))

		y,a,b = max(Y)
		
		sel = {a:None,b:None}
		
		LV[tuple(sorted(list(np.append(nod[a],nod[b]))))] = (nod[a],nod[b])

		nod = [nod[i] for i in range(len(nod)) if not sel.has_key(i)] + [np.append(nod[a],nod[b])]
	
	return LV
	
def Boot(X,Nt=1000):
	Rb = []
	for i in range(Nt):
		sel = np.random.choice(range(X.shape[1]),replace=True,size=X.shape[1])
		Xb = X[:,sel]
		Rb.append(np.corrcoef(Xb))
	return Rb

def Collect(x,Rb,LV):
    d = [Rb[i][LV[srt(x)][0]][:,LV[srt(x)][1]].mean() for i in range(len(Rb))]
    return np.array(d)


def Compare_Null(LV,Rb):
	C = {k:Collect(k,Rb,LV) for k in LV}

	P = []
	KEY = C.keys()
	for i in range(len(KEY)):
		for j in range(i):
			if len(np.intersect1d(KEY[i],KEY[j]))==0: continue
			if not (KEY[i] in map(tuple,map(sorted,LV[KEY[j]])) or KEY[j] in map(tuple,map(sorted,LV[KEY[i]]))): continue
			p= getPvalue(C[KEY[i]],C[KEY[j]])
			
			P.append((p,i,j))

	P = np.array(sorted(P))
	return KEY,P

def Find_Component(LV,Rb,alpha=0.05):
	
	KEY,P = Compare_Null(LV,Rb):
	indx = np.where(P[:,0]<alpha*np.arange(1,P.shape[0]+1)/P.shape[0])[0]

	if len(indx)==0:
		indx = 0
	else:
		indx = indx[-1]+1


	ed = list(map(tuple,P[indx:,1:].astype(int)))
	g = ig.Graph(n=len(C),edges=ed)
	g.vs["name"] = KEY

	comp = g.components()

	return comp


def HValidate(LV,Rb,alpha=0.05):
	comp = Find_Component(LV,RB,alpha)

	comp = comp.membership
	L = [[] for i in range(len(set(comp)))]
	for i,c in enumerate(comp):
		L[c].extend(KEY[i])
	L = [list(set(l)) for l in L]

	return L
	
def HclustVal(X,alpha=0.05,Nt=1000):
	R = np.corrcoef(X)
	
	LV = AverageLink(R)
	
	Rb = Boot(X,Nt)
	
	L = HValidate(LV,RB,alpha)
	
	return L,ToMemb(L,X.shape[0])
	
	
