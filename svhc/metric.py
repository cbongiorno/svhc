import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import scipy.special as ss
from collections import defaultdict,OrderedDict
import matplotlib.pyplot as plt
#from MainLib import ToMemb
import sys
import pandas as pd

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

def InterCorr(R,C,g):
    
    N = R.shape[0]
    A = np.zeros((N,N))
    
    v = np.zeros(N)
    v[np.array(g)] = 1
    part = np.outer(v,v)
    A = A+part
    
    for cl in C:
        cl = np.array(cl)
        v = np.zeros(N)
        v[cl] = 1

        part = np.outer(v,v)
        A = A-part
    np.fill_diagonal(A,0)
    return np.mean(R[np.where(A)])

def Find_Gen(L):
    gen = {}
    for i in range(len(L)):
        for j in range(i+1,len(L)):
            if len(set(L[i]) & set(L[j]))==len(L[i]):
                gen[L[i]] = L[j]
                break
    return gen

def CreateAvMat(X,L):
    N = X.shape[0]
    R = np.corrcoef(X)
    La = [tuple(sorted(l)) for l in L]+[(i,) for i in range(N)]

    La = sorted(La,key=len)

    gen = Find_Gen(La)

    son = defaultdict(list)
    for s,g in gen.items():
        son[g].append(s)

    La =sorted(La,key=len,reverse=True)


    Rm = np.zeros((N,N))

    for i in range(len(La)):
        if not son.has_key(La[i]): continue
        #if len(La[i])==1: continue
        p = InterCorr(R,son[La[i]],La[i])
        #print p,i
        v = np.zeros(N)
        v[np.array(La[i])] = 1
        v = np.outer(v,v)
        Rm[np.where(v)] = p

    return Rm

def ComputeDistance(Rx,Rb):
    return np.mean(map(lambda x:abs(x-Rx)[np.triu_indices(Rx.shape[0],1)].mean(),Rb))

def Performance(Rx,R,Rb):
    mx = ComputeDistance(R,Rb)

    return (ComputeDistance(Rx,Rb)-mx)/(ComputeDistance(np.identity(R.shape[0]),Rb)-mx)

def _ARIparts(m1,m2):
    M1 = defaultdict(list)
    M2 = defaultdict(list)
    N = len(m1)
    for i in xrange(len(m1)):
        M1[m1[i]].append(i)
        M2[m2[i]].append(i)


    A = sum([ss.binom(len(set(M1[i]).intersection(M2[j])),2) for i in M1 for j in M2])
    B = sum([ss.binom(sum([len(set(M1[i]).intersection(M2[j])) for j in M2]),2) for i in M1])
    C = sum([ss.binom(sum([len(set(M1[i]).intersection(M2[j])) for i in M1]),2) for j in M2])
    return A,B,C

def HierARI(M1,M2):
    M1 = M1[:]
    M2 = M2[:]
    
    if len(M1)<len(M2):
        for i in range(len(M2)-len(M1)):
            M1.append(M1[-1])
    if len(M2)<len(M1):
        for i in range(len(M1)-len(M2)):
            M2.append(M2[-1])
    M1 = np.array(M1)
    M2 = np.array(M2)
    M,N = M1.shape
    P = [_ARIparts(M1[i],M2[i]) for i in xrange(M)]
    A,B,C = map(np.array,zip(*P))
    
    return (A - B*C/ss.binom(N,2)).mean()/(0.5*(B+C)-B*C/ss.binom(N,2)).mean()
   
def _awiPart(A,B):
    
    mb = defaultdict(list)
    for i in xrange(len(A)):
        mb[A[i]].append(i)
    mb = [mb[e] for e in mb if len(mb[e])>1]

    mr = defaultdict(list)
    for i in xrange(len(B)):
        mr[B[i]].append(i)
    mr = [mr[e] for e in mr if len(mr[e])>1]
    Bp = sum(map(lambda x:len(x)*(len(x)-1)/2,mr))

    N = len(A)*(len(A)-1)/2

    n = 0
    m = 0
    for x in mb:
        for i in xrange(len(x)):
            for j in range(i+1,len(x)):
                n+=1
                if B[x[i]]==B[x[j]]:
                    m+=1

    Ap = float(n)
    num = m
    den = Ap
    exp = (Bp*Ap)/N
    return num,den,exp
    
def HierAWI(M1,M2):
    M1 = M1[:]
    M2 = M2[:]
    if len(M1)<len(M2):
        for i in range(len(M2)-len(M1)):
            M1.append(M1[-1])
    if len(M2)<len(M1):
        for i in range(len(M1)-len(M2)):
            M2.append(M2[-1])
    M1 = np.array(M1)
    M2 = np.array(M2)
    
    P = np.array([_awiPart(M1[i],M2[i]) for i in xrange(M1.shape[0])])
    
    num,den,exp = P.sum(axis=0)
    if (den - exp) == 0:
		return None
		
    return (num-exp)/float(den-exp)

    
def StandardDendrogram(H):
	h = sorted(zip(*H))
	h = map(list,h)

	for j in range(len(h[0])-1):
		k = OrderedDict.fromkeys(zip(*h)[j]).keys()
		k = dict(zip(k,range(len(k))))

		for i in range(len(h)):
			h[i][j] = k[h[i][j]]

	j = len(h[0])-1

	k = OrderedDict.fromkeys(zip(*h)[j]).keys()
	k = dict(zip(k,range(len(k))))

	h = sorted(h,key=lambda x:x[-1])

	for i in range(len(h)):
			h[i][j] = k[h[i][j]]
			
	h = map(np.array,zip(*h))
	return h

def OrderArray(H):
    return np.array( zip(*sorted(zip(H[-1],range(len(H[0])))))[1] )

'''
def CompressHC(nodes,N):
	return StandardDendrogram(ToMemb(nodes,N))
'''
	
if __name__=='__main__':
	
	partX,partRef = sys.argv[1:]
	
	partX = map(np.array,np.array(pd.read_csv(partX,sep='\t',header=None)))
	partRef = map(np.array,np.array(pd.read_csv(partRef,sep='\t',header=None)))
	
	hari,hawi = HierARI(partX,partRef),HierAWI(partX,partRef)
	
	print "HARI",hari
	print "HAWI",hawi
	
	with open('metric_out.dat','w') as fw:
		fw.write("%f\tHARI\n%f\tHAWI\n"%(hari,hawi))
	
	
	
