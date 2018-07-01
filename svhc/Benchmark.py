import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import scipy.stats as st
import sys
import pandas as pd

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

	
	
