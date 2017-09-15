import sys
import hclustval
import numpy as np
import pandas as pd


if __name__=='__main__':
	inputfile,T,noise = sys.argv[1:]
	T,noise = int(T),float(noise)
	
	P = np.array(pd.read_csv(inputfile,sep='\t',header=None))
	
	X,comm = hclustval.CREATE_DATA(P,T=T,noise=noise)
	
	pd.DataFrame(X).to_csv('dataSeries_benchmark.dat',sep='\t',index=False,header=False)
	with open('HyperNodeList_reference.dat','w') as fw:
		fw.write("\n".join([",".join(map(str,comm[i])) for i in range(len(comm))]))
	
	pd.DataFrame(hclustval.CompressHC(comm,X.shape[0])).to_csv('CompressDendrogram_reference.dat',sep='\t',index=False,header=False)
	
