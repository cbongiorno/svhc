import matplotlib.pyplot as plt
import matplotlib.cm as cmap
from metric import StandardDendrogram,OrderArray
from MainLib import ToMemb
import matplotlib.patches as patches


def PlotDendro(R,nodes,xylabel='Objects',cblabel='pearson',file_w=None,lw=0.5,color='red'):

	ordn = OrderArray(StandardDendrogram(ToMemb(nodes,R.shape[0])))
	rev_ornd = zip(*sorted(zip(ordn,range(len(ordn)))))[1]

	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111, aspect='equal')

	for l in nodes:
		l = sorted([rev_ornd[a] for a in l])

		ax2.add_patch(
			patches.Rectangle(
				(min(l)-0.5, min(l)-0.6),
				len(l),
				len(l),
				fill=False,
				lw=lw,
				color=color,	
				)
			)
		

	plt.imshow(R[ordn][:,ordn],vmin=-1,vmax=1,cmap=cmap.RdBu_r,interpolation='nearest')
	plt.ylabel(xylabel,fontsize=16)
	plt.xlabel(xylabel,fontsize=16)
	plt.xticks(fontsize=11)
	plt.yticks(fontsize=11)
	cb = plt.colorbar()
	cb.set_label(cblabel,fontsize=14)
	plt.tight_layout()
	if file_w!=None:
		plt.savefig(file_w)
	plt.show()
