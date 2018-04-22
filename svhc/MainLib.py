import numpy as np
import igraph as ig
import sys
import pandas as pd
import hclustval


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

def SingleLink(R):
	N = R.shape[0]
	LV = {}

	nod =  [np.array([i]) for i in xrange(N)]

	for lv in xrange(N-1):
		Y = []
		for i in xrange(len(nod)):
			for j in xrange(i):
				Y.append((R[nod[i]][:,nod[j]].max(),i,j))

		y,a,b = max(Y)
		
		sel = {a:None,b:None}
		
		LV[tuple(sorted(list(np.append(nod[a],nod[b]))))] = (nod[a],nod[b])

		nod = [nod[i] for i in range(len(nod)) if not sel.has_key(i)] + [np.append(nod[a],nod[b])]
	
	return LV
	
def CompleteLink(R):
	N = R.shape[0]
	LV = {}

	nod =  [np.array([i]) for i in xrange(N)]

	for lv in xrange(N-1):
		Y = []
		for i in xrange(len(nod)):
			for j in xrange(i):
				Y.append((R[nod[i]][:,nod[j]].min(),i,j))

		y,a,b = max(Y)
		
		sel = {a:None,b:None}
		
		LV[tuple(sorted(list(np.append(nod[a],nod[b]))))] = (nod[a],nod[b])

		nod = [nod[i] for i in range(len(nod)) if not sel.has_key(i)] + [np.append(nod[a],nod[b])]
	
	return LV

def pairBoot(Rb):    
    Nt,N = len(Rb),Rb[0].shape[0]
    
    for i in xrange(N):
        for j in xrange(i):    
            ordn = np.random.choice(Nt,size=Nt,replace=False)
            cmd = [Rb[x][i,j] for x in xrange(Nt)]
            for x in xrange(Nt):
                Rb[x][i,j] = cmd[x]
                Rb[x][j,i] = cmd[x]
    return Rb	
	
def Boot(X,Nt=1000,boottype='standard',pattern=None,M=None):
	
	if pattern==None:
		Rb = []
		for i in xrange(Nt):
			sel = np.random.choice(range(X.shape[1]),replace=True,size=X.shape[1])
			Xb = X[:,sel]
			Rb.append(np.corrcoef(Xb))

		if boottype=='pair':
			Rb = pairBoot(Rb)
		return Rb
	else:
		Rb = []
		for _ in xrange(Nt):
			Xb,cmm = hclustval.Benchmark.CREATE_DATA(pattern,T=M)
			Rb.append(np.corrcoef(Xb))
			
		return Rb

def Collect(x,Rb,LV,method):
	
	if method=='average':
		d = [Rb[i][LV[srt(x)][0]][:,LV[srt(x)][1]].mean() for i in range(len(Rb))]
	elif method=='single':
		d = [Rb[i][LV[srt(x)][0]][:,LV[srt(x)][1]].max() for i in range(len(Rb))]
	elif method=='complete':
		d = [Rb[i][LV[srt(x)][0]][:,LV[srt(x)][1]].min() for i in range(len(Rb))]
		
	return np.array(d)


def Compare_Null(LV,Rb,method):
	C = {k:Collect(k,Rb,LV,method) for k in LV}

	P = []
	KEY = C.keys()
	for i in xrange(len(KEY)):
		for j in xrange(i):
			if len(np.intersect1d(KEY[i],KEY[j]))==0: continue
			if not (KEY[i] in map(tuple,map(sorted,LV[KEY[j]])) or KEY[j] in map(tuple,map(sorted,LV[KEY[i]]))): continue
			p= getPvalue(C[KEY[i]],C[KEY[j]])
			
			P.append((p,i,j))

	P = np.array(sorted(P))
	return KEY,P

def Find_Component(LV,Rb,method,alpha=0.05):
	
	KEY,P = Compare_Null(LV,Rb,method)
	indx = np.where(P[:,0]<alpha*np.arange(1,P.shape[0]+1)/P.shape[0])[0]

	if len(indx)==0:
		indx = 0
	else:
		indx = indx[-1]+1


	ed = list(map(tuple,P[indx:,1:].astype(int)))
	g = ig.Graph(n=len(KEY),edges=ed)
	g.vs["name"] = KEY

	comp = g.components()

	return comp


def HValidate(LV,Rb,method,alpha=0.05):
	comp = Find_Component(LV,Rb,method,alpha)

	KEY = comp.graph.vs["name"]
	comp = comp.membership
	L = [[] for i in range(len(set(comp)))]
	for i,c in enumerate(comp):
		L[c].extend(KEY[i])
	L = [list(set(l)) for l in L]

	return L
	
def HclustVal(X,alpha=0.05,Nt=1000,method='average',boottype='standard',Rb=None,pattern=None,M=None):
	
	if pattern!=None and M==None:
		print "If you specify the pattern, you must specify the dataseries lenght (M) too"
		sys.exit(0)
	
	
	R = np.corrcoef(X)
	
	if method!='all':
		if method=='average':
			LV = AverageLink(R)
		elif method=='single':
			LV = SingleLink(R)
		elif method=='complete':
			LV = CompleteLink(R)
		
		if Rb==None:
			Rb = Boot(X,Nt,boottype,pattern,M)
		
		L = HValidate(LV,Rb,method,alpha)
		
		return L,hclustval.StandardDendrogram(ToMemb(L,X.shape[0])),LV
	
	else:
		res = {}
		
		res['single'] = SingleLink(R)
		if Rb==None:
			Rb = Boot(X,Nt,boottype,pattern,M)
		res['singleVal'] = HValidate(res['single'],Rb,'single',alpha)
		#res['singValCompress'] = hclustval.StandardDendrogram(ToMemb(res['single'],X.shape[0]))
		
		res['average'] = AverageLink(R)
		if Rb==None:
			Rb = Boot(X,Nt,boottype,pattern,M)
		res['averageVal'] = HValidate(res['average'],Rb,'average',alpha)
		#res['averageValCompress'] = hclustval.StandardDendrogram(ToMemb(res['average'],X.shape[0]))
		
		res['complete'] = CompleteLink(R)
		if Rb==None:
			Rb = Boot(X,Nt,boottype,pattern,M)
		res['completeVal'] = HValidate(res['complete'],Rb,'complete',alpha)
		#res['completeValCompress'] = hclustval.StandardDendrogram(ToMemb(res['complete'],X.shape[0]))
		
		return res
		
		
	
if __name__=='__main__':
	
	inputfile,alpha,Nt,method = sys.argv[1:]
	
	alpha,Nt = float(alpha),int(Nt)
	
	X = np.array(pd.read_csv(inputfile,sep='\t',header=None))
	
	comm,H = HclustVal(X,alpha,Nt,method)
	
	with open('HyperNodeList_out.dat','w') as fw:
		fw.write("\n".join([",".join(map(str,comm[i])) for i in range(len(comm))]))

	pd.DataFrame(hclustval.CompressHC(comm,X.shape[0])).to_csv('CompressDendrogram_out.dat',sep='\t',index=False,header=False)


	
	
	
