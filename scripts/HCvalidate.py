import sys
import hclustval
import pandas as pd
import numpy as np

if __name__=='__main__':
	inputfile,alpha,Nt,method = sys.argv[1:]
	
	alpha,Nt = float(alpha),int(Nt)
	
	X = np.array(pd.read_csv(inputfile,sep='\t',header=None))
	
	comm,H = hclustval.HclustVal(X,alpha,Nt,method)
	H = hclustval.StandardDendrogram(H)
	
	with open('HyperNodeList_out.dat','w') as fw:
		fw.write("\n".join([",".join(map(str,comm[i])) for i in range(len(comm))]))

	pd.DataFrame(hclustval.CompressHC(comm,X.shape[0])).to_csv('CompressDendrogram_out.dat',sep='\t',index=False,header=False)

