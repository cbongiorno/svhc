import scipy.stats as st
import numpy as np
import pandas as pd
import sys

def GarlaschelliBench(Q,N,M,mu,ni):
    a = st.norm.rvs(0,1,size=M)
    b = st.norm.rvs(0,1,size=(N,M))
    g = st.norm.rvs(0,1,size=(len(Q),M))

    a = np.array([a for i in xrange(N)])

    g = np.array([g[j] for j in xrange(len(Q)) for i in xrange(Q[j])])

    X = mu*a + ni*b + g
    return st.zscore(X,axis=1)
    
if __name__=='__main__':
	
	inputfile,mu,ni,M = sys.argv[1:]
	
	mu,ni,M = float(mu),float(ni),int(M)
	
	Q = np.array(pd.read_csv(inputfile,header=None,sep=','))[0]
	
	X = GarlaschelliBench(Q,Q.sum(),M,mu,ni)
	
	pd.DataFrame(X).to_csv('dataSeries_benchmark.dat',sep='\t',index=False,header=False)
	
	H = [[0]*Q.sum(),[j for j in xrange(len(Q)) for i in xrange(Q[j])],range(Q.sum())]	
	pd.DataFrame(H).to_csv('CompressDendrogram_reference.dat',sep='\t',index=False,header=False)
	
	com = []
	h = 0
	for i in xrange(Q.shape[0]):
		a = []
		for j in xrange(Q[i]):
			a.append(h)
			h+=1
		com.append(a)
			

	com.append(range(Q.sum()))
	com = sorted(com,key=lambda x:len(x),reverse=True)
	
	with open('HyperNodeList_reference.dat','w') as fw:
		fw.write("\n".join([",".join(map(str,com[i])) for i in range(len(com))]))
