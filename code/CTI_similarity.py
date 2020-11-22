# calculate pairwise chemical-chemical CTI similarity given a list of chemicals

import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd 

folder = '../data/'
for tag in ['Antiviral/','Anticytokine/']: 
	chem_file = folder + tag + 'top_compounds.txt'
	t_file = folder + tag + 'targets_for_compounds.txt'
	dt_file = folder + tag + 'interactions_for_compounds.txt'
	simi_out = folder + tag + 'compounds_with_CTI_similarity.txt'

	chem_df = pd.read_csv(chem_file, sep = '\t', usecols = ['Index'])	
	cids = chem_df['Index'].tolist()
	cid2ind = {}
	for i, cid in enumerate(cids):
		cid2ind[cid] = i

	pro_df = pd.read_csv(t_file, sep = '\t', usecols = ['Uniprot_ID'])
	unips = pro_df['Uniprot_ID'].tolist()
	unip2ind = {}
	for i, unip in enumerate(unips):
		unip2ind[unip] = i

	dt_df = pd.read_csv(dt_file, sep = '\t', usecols = ['Index', 'Uniprot_ID', 'Confidence_Score'])
	
	dt_matrix = np.zeros(shape=(len(cids), len(unips)))
	for _, row in dt_df.iterrows():
		cid, unip, score = row['Index'], row['Uniprot_ID'], row['Confidence_Score']
		i, j = cid2ind[cid], unip2ind[unip]
		dt_matrix[i][j] = score

	s = cosine_similarity(dt_matrix,dt_matrix)

	out = open(simi_out, 'w')
	out.write('Index_1\tIndex_2\tSimilarity_CTI\n')
	n = len(cids)
	for i in range(n):
		for j in range(i + 1, n):
			s1 = str(s[i][j])
			outline = [str(cids[i]), str(cids[j]), s1]
			out.write('\t'.join(outline) + '\n')
	out.close()