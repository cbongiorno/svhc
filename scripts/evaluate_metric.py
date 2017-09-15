import sys
import hclustval
import numpy as np
import pandas as pd
	
if __name__=='__main__':
	
	partX,partRef = sys.argv[1:]
	
	partX = map(np.array,np.array(pd.read_csv(partX,sep='\t',header=None)))
	partRef = map(np.array,np.array(pd.read_csv(partRef,sep='\t',header=None)))
	
	hari,hawi = hclustval.HierARI(partX,partRef),hclustval.HierAWI(partX,partRef)
	
	print "HARI",hari
	print "HAWI",hawi
	
	with open('metric_out.dat','w') as fw:
		fw.write("%f\tHARI\n%f\tHAWI\n"%(hari,hawi))
	
