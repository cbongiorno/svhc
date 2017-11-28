import matplotlib.pyplot as plt
import matplotlib.cm as cmap
from metric import StandardDendrogram,OrderArray
from MainLib import ToMemb
import matplotlib.patches as patches
import numpy as np
import seaborn as sns

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

#L,H,LV = hclustval.HclustVal(x)

from hclustval.metric import OrderArray,StandardDendrogram,ToMemb



def DendroAndCorrDist(x,L,LV,method='average',xylabel='Objects',file_w=None,color_style='Set1'):
    
	ordn = ordn = OrderArray(StandardDendrogram(ToMemb(LV.keys(),x.shape[0])))
	x = x[ordn]

	rev_ornd = dict(sorted(zip(ordn,range(len(ordn)))))

	LV = {tuple(sorted([rev_ornd[k] for k in K])):(np.array(sorted([rev_ornd[a] for a in LV[K][0]])),
												   np.array(sorted([rev_ornd[a] for a in LV[K][1]]))) for K in LV}

	L = [tuple(sorted([rev_ornd[k] for k in K])) for K in L]
	R = 1.-np.corrcoef(x)

	L = sorted(L,key=len)

	N = R.shape[0]
	K = LV.keys()

	import matplotlib.patches as patches
	import colorsys
	nc  = len(L)
	import matplotlib as mpl
	import matplotlib.cm as cm

	cnmx = sns.color_palette(color_style, nc-1)

	from collections import defaultdict
	col = defaultdict(int)
	for k in K:
		for j in range(len(L)):
			if set(k) <= set(L[j]):

				if tuple(L[j])==tuple(range(N)):
					col[tuple(k)] = 'black'
					break

				col[tuple(k)] = cnmx[j-1]
				break


	from matplotlib import gridspec

	fig = plt.figure(figsize=(5,9))
	gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.34]) 

	# the fisrt subplot
	axarr = [0,0]
	axarr[0] = plt.subplot(gs[0])
	axarr[1] = plt.subplot(gs[1], sharex = axarr[0])


	pos = defaultdict(tuple)

	for i in range(N):
		pos[(i,)] = (float(i),1)

	K = sorted(LV.keys(),key=len)

	hm = []
	for k in K:
		d = tuple(sorted(LV[k][0])),tuple(sorted(LV[k][1]))

		c = tuple(sorted(k))
		#if not col.has_key(c):
		#    break
		pp = pos[d[0]][0],pos[d[1]][0]
		h0 = pos[d[0]][1],pos[d[1]][1]

		if method=='average':
			h = R[LV[k][0]][:,LV[k][1]].mean()
		elif method=='complete':
			h = R[LV[k][0]][:,LV[k][1]].min()
		elif method=='single':
			h = R[LV[k][0]][:,LV[k][1]].max()

		axarr[0].plot([pp[0],pp[0]],[h0[0],h],'-',color='k',alpha=0.5)
		axarr[0].plot([pp[1],pp[1]],[h0[1],h],'-',color='k',alpha=0.5)
		axarr[0].plot([pp[0],pp[1]],[h,h],'-',color='k',lw=1.5)
		axarr[0].plot([pp[0],pp[1]],[h,h],'-',color=col[c],lw=1.4,zorder=1000)

		hm.append(h)

		#plt.plot((pp[0]+pp[1])/2,h,'o',color=col[c],ms=4,zorder=1000)

		pos[c] = (pp[0]+pp[1])/2,h
		#print "ciao"

	axarr[0].set_yticks(np.arange(0,1.2,0.2))
	axarr[0].set_ylim([0,max(hm)+0.05])

	axarr[0].yaxis.tick_right()
	axarr[0].tick_params(labelsize=14)
	axarr[0].set_ylabel('Metric (%s)'%method,fontsize=18)
	axarr[0].yaxis.set_label_position("right")

	ims = axarr[1].imshow(R,vmin=-1,vmax=1,cmap=cm.RdBu_r,interpolation='nearest')


	cbar_ax = fig.add_axes([0.95, 0.133, 0.03, 0.4])
	cbr = fig.colorbar(ims, cax=cbar_ax)

	cbr.set_label('correlation',fontsize=18)
	axarr[1].set_xlabel(xylabel,fontsize=18)
	axarr[1].set_ylabel(xylabel,fontsize=18)
	axarr[1].tick_params(labelsize=14)

	axarr[1].set_ylim([0,N])

	for l in L:
		if l == tuple(range(N)): continue

		axarr[1].add_patch(
			patches.Rectangle(
				(min(l)-0.5, min(l)-0.5),
				len(l),
				len(l),
				fill=False,
				lw=0.8,
				color='green',	
				)
			)

	plt.setp(axarr[0].get_xticklabels(), visible=False)
	plt.subplots_adjust(hspace=.0)

	if file_w!=None:
		plt.savefig('%s_%s.pdf'%(file_w,method), bbox_extra_artists=(cbar_ax,), bbox_inches='tight')

	plt.show()
	return



def DendroAndCorr(x,L,LV,method='average',xylabel='Objects',file_w=None,color_style='Set1'):
    
	ordn = ordn = OrderArray(StandardDendrogram(ToMemb(LV.keys(),x.shape[0])))
	x = x[ordn]

	rev_ornd = dict(sorted(zip(ordn,range(len(ordn)))))

	LV = {tuple(sorted([rev_ornd[k] for k in K])):(np.array(sorted([rev_ornd[a] for a in LV[K][0]])),
												   np.array(sorted([rev_ornd[a] for a in LV[K][1]]))) for K in LV}

	L = [tuple(sorted([rev_ornd[k] for k in K])) for K in L]
	R = np.corrcoef(x)

	L = sorted(L,key=len)

	N = R.shape[0]
	K = LV.keys()

	import matplotlib.patches as patches
	import colorsys
	nc  = len(L)
	import matplotlib as mpl
	import matplotlib.cm as cm

	cnmx = sns.color_palette(color_style, nc-1)

	from collections import defaultdict
	col = defaultdict(int)
	for k in K:
		for j in range(len(L)):
			if set(k) <= set(L[j]):

				if tuple(L[j])==tuple(range(N)):
					col[tuple(k)] = 'black'
					break

				col[tuple(k)] = cnmx[j-1]
				break


	from matplotlib import gridspec

	fig = plt.figure(figsize=(5,9))
	gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.34]) 

	# the fisrt subplot
	axarr = [0,0]
	axarr[0] = plt.subplot(gs[0])
	axarr[1] = plt.subplot(gs[1], sharex = axarr[0])


	pos = defaultdict(tuple)

	for i in range(N):
		pos[(i,)] = (float(i),1)

	K = sorted(LV.keys(),key=len)

	hm = []
	for k in K:
		d = tuple(sorted(LV[k][0])),tuple(sorted(LV[k][1]))

		c = tuple(sorted(k))
		#if not col.has_key(c):
		#    break
		pp = pos[d[0]][0],pos[d[1]][0]
		h0 = pos[d[0]][1],pos[d[1]][1]

		if method=='average':
			h = R[LV[k][0]][:,LV[k][1]].mean()
		elif method=='complete':
			h = R[LV[k][0]][:,LV[k][1]].min()
		elif method=='single':
			h = R[LV[k][0]][:,LV[k][1]].max()

		axarr[0].plot([pp[0],pp[0]],[h0[0],h],'-',color='k',alpha=0.5)
		axarr[0].plot([pp[1],pp[1]],[h0[1],h],'-',color='k',alpha=0.5)
		axarr[0].plot([pp[0],pp[1]],[h,h],'-',color='k',lw=1.5)
		axarr[0].plot([pp[0],pp[1]],[h,h],'-',color=col[c],lw=1.4,zorder=1000)

		hm.append(h)

		#plt.plot((pp[0]+pp[1])/2,h,'o',color=col[c],ms=4,zorder=1000)

		pos[c] = (pp[0]+pp[1])/2,h
		#print "ciao"

	axarr[0].set_yticks(np.arange(0,1.2,0.2))
	axarr[0].set_ylim([1,min(hm)-0.05])

	axarr[0].yaxis.tick_right()
	axarr[0].tick_params(labelsize=14)
	axarr[0].set_ylabel('Metric (%s)'%method,fontsize=18)
	axarr[0].yaxis.set_label_position("right")

	ims = axarr[1].imshow(R,vmin=-1,vmax=1,cmap=cm.RdBu_r,interpolation='nearest')


	cbar_ax = fig.add_axes([0.95, 0.133, 0.03, 0.4])
	cbr = fig.colorbar(ims, cax=cbar_ax)

	cbr.set_label('correlation',fontsize=18)
	axarr[1].set_xlabel(xylabel,fontsize=18)
	axarr[1].set_ylabel(xylabel,fontsize=18)
	axarr[1].tick_params(labelsize=14)

	axarr[1].set_ylim([0,N])

	for l in L:
		if l == tuple(range(N)): continue

		axarr[1].add_patch(
			patches.Rectangle(
				(min(l)-0.5, min(l)-0.5),
				len(l),
				len(l),
				fill=False,
				lw=0.8,
				color='green',	
				)
			)

	plt.setp(axarr[0].get_xticklabels(), visible=False)
	plt.subplots_adjust(hspace=.0)

	if file_w!=None:
		plt.savefig('%s_%s.pdf'%(file_w,method), bbox_extra_artists=(cbar_ax,), bbox_inches='tight')

	plt.show()
	return
