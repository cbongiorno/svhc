import numpy as np
import scipy.stats as st
import sys
import pandas as pd
import hclustval

def _gs(X, row_vecs=True, norm = True):
    if not row_vecs:
        X = X.T
    Y = X[0:1,:].copy()
    for i in range(1, X.shape[0]):
        proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, X[i,:] - proj.sum(0)))
    if norm:
        Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
    if row_vecs:
        return Y
    else:
        return Y.T

def CREATE_DATA(P,T=1000,noise=0.,com_ortho=True):
    N = P.shape[0]

    U = np.diag(np.sqrt(1-(P**2).sum(axis=1)))

    B = np.hstack((P,U))

    A = st.norm.rvs(0,1,size=(B.shape[1],T))
    if com_ortho==True:
		A[:P.shape[1]] = _gs(A[:P.shape[1]])
    A =  A/np.linalg.norm(A,axis=1)[np.newaxis].T

    X = np.dot(B,A)
    
    com = map(lambda x:np.where(x)[0],(P>0).T)
    X = st.zscore(X,axis=1)
    
    X = (1.-noise)*X + noise*st.norm.rvs(0,1,size=X.shape)
    
    return X,com

if __name__=='__main__':
	
	inputfile,T,noise = sys.argv[1:]
	T,noise = int(T),float(noise)
	
	P = np.array(pd.read_csv(inputfile,sep='\t',header=None))
	
	X,comm = CREATE_DATA(P,T=T,noise=noise)
	
	pd.DataFrame(X).to_csv('dataSeries_benchmark.dat',sep='\t',index=False,header=False)
	with open('HyperNodeList_reference.dat','w') as fw:
		fw.write("\n".join([",".join(map(str,comm[i])) for i in range(len(comm))]))
	
	pd.DataFrame(hclustval.CompressHC(comm,X.shape[0])).to_csv('CompressDendrogram_reference.dat',sep='\t',index=False,header=False)
	
	
