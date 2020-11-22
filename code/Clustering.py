#cluster chemicals based on their similarities

import pandas as pd
import numpy as np 
import matplotlib as mpl
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
from scipy.stats import hypergeom
import matplotlib.pyplot as plt

def get_sorted_ind(x):
	"""
	Takes a item-item distance matrix, cluster it based on the linkage distance, return the sorted row or column index
	"""
	d = dist.pdist(x)
	D = dist.squareform(d)
	Y = sch.linkage(D, method='average', metric='cosine') 
	Z = sch.dendrogram(Y)
	idx = Z['leaves'] 
	return idx	

def heatmap(x, xlabel, title, filename):
	"""
	This below code is based in large part on the protype methods:
	http://old.nabble.com/How-to-plot-heatmap-with-matplotlib--td32534593.html
	http://stackoverflow.com/questions/7664826/how-to-get-flat-clustering-corresponding-to-color-clusters-in-the-dendrogram-cre
	"""	
	cmap=plt.cm.bwr
	#cmap=plt.get_cmap('coolwarm')
	default_window_hight = 10
	default_window_width = 10
	fig = plt.figure(figsize=(default_window_width,default_window_hight))
	axm = fig.add_axes([0.15, 0.2, 0.65, 0.65])
	im = axm.matshow(x, aspect='auto', origin='lower', cmap=cmap)
	m = len(x)
	if m <=20:
		step =1
	else:
		step = round((m-1)/20)+1
		#step = 6
	axm.xaxis.set_ticks_position('bottom')
	axm.set_xticks(np.arange(0,m,step),minor=False)
	axm.set_xticklabels(np.arange(1,m+1,step),minor=False)
	axm.set_yticks(np.arange(0,m,step),minor=False)
	axm.set_yticklabels(np.arange(1,m+1,step),minor=False)
	axm.set_xlabel('\n'+xlabel,fontsize= 18,fontweight='bold')
	axm.set_ylabel(xlabel+'\n',fontsize= 18,fontweight='bold')
	axm.set_title(title,fontsize= 24,fontweight='bold')

	axcb = fig.add_axes([0.85, 0.2, 0.02, 0.65], frame_on=False)  # axes for colorbar
	cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, orientation='vertical')
	font_size = 12 # Adjust as appropriate.
	cb.ax.tick_params(labelsize=font_size)
	#axcb.set_title("Smilarity\n",fontsize= 16)

	plt.rcParams['font.size'] = 16

	plt.savefig(filename, dpi=800)
	plt.close()

def write_matrix(x, row_header, column_header, filename):
	out = open(filename,'w')
	out.write('\t'+'\t'.join(column_header)+'\n')
	i = 0
	for row in x:
		out.write(row_header[i]+'\t'+'\t'.join(map(str,row))+'\n')
		i+=1
	out.close()

def get_clusters(simi_file, cid2ind, cids, x_label, title, simi_plot, out_matrix):
	n = len(cid2ind)
	simi_x = np.zeros(shape=(n, n))
	for i in range(n):
		simi_x[i][i] = 1.0

	with open(simi_file, 'r') as f:
		next(f)
		for line in f:
			line = line.strip().split('\t')
			if line[0] not in cid2ind or line[1] not in cid2ind:
				continue
			i, j, s = cid2ind[line[0]], cid2ind[line[1]], float(line[2])
			simi_x[i][j] = s
			simi_x[j][i] = s

	idx = get_sorted_ind(simi_x)
	sorted_cids = [cids[k] for k in idx]
	sorted_X = simi_x[idx]
	sorted_X = sorted_X[:, idx]
	heatmap(sorted_X, x_label, title, simi_plot)
	write_matrix(sorted_X, sorted_cids, sorted_cids, out_matrix)

	cid2new_order = {sorted_cids[i]: i for i in range(len(cids))}

	return cid2new_order

def my_main():
	folder = '../data/'
	for tag in ['Antiviral/clustering/','Anticytokine/clustering/']: 
		chem_file = folder + tag + 'top_compounds.txt'
		chem_out = folder + tag + 'top_compounds_with_clusters.txt'
		x_label = 'Prioritized compounds'

		chem_df = pd.read_csv(chem_file, sep = '\t')
		cids = chem_df['Index'].tolist()
		cids = [str(x) for x in cids]
		cid2ind = {}
		for i, cid in enumerate(cids):
			cid2ind[cid] = i

		simi_file = folder + tag + 'compounds_with_CTI_similarity.txt'
		simi_plot = folder + tag + 'compounds_with_CTI_similarity.png'
		out_matrix = folder + tag + 'compounds_with_CTI_similarity_matrix.txt'
		title = 'Interaction-pattern-based clustering' 

		cid2new_order = get_clusters(simi_file, cid2ind, cids, x_label, title, simi_plot, out_matrix)
		chem_df['order_CTI'] = chem_df.apply(lambda x: cid2new_order[str(x['Index'])], axis = 1)
	chem_df.fillna('NA', inplace = True)
	chem_df.to_csv(chem_out, index = False, sep = '\t')

if __name__ == '__main__':
	my_main()
