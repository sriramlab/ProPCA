#!/bin/env python3

import sys
import numpy as np
import pandas as pd
from pandas_plink import read_plink
import argparse

if __name__=='__main__':

	parser = argparse.ArgumentParser(description="Calculate chi-squared selection statistics based on principal components from Galinsky et al. 2016")
	parser.add_argument("outfile",help="path to output file name")
	parser.add_argument("plink",help="path to PLINK prefix")
	parser.add_argument("eval",help="path to eigenvalue file")
	parser.add_argument("proj",help="path to projections file")
	parser.add_argument("-v","--verbose",help="verbose mode (default: TRUE)",action="store_false")
	parser.add_argument("-m","--missing",help="missing mode (default: FALSE)",action="store_true")
	parser.add_argument("-c","--chunk",help="chunk size (default: 64)",type=int,default=64)

	args = parser.parse_args()
	
	outfile = args.outfile
	filename = args.plink
	eigenvec_file = args.proj
	eigenvals_file = args.eval	
	verbose = args.verbose
	chunk_size = args.chunk
	missing = args.missing

	evecs = np.loadtxt(eigenvec_file,dtype=np.float64)
	evals = np.loadtxt(eigenvals_file,dtype=np.float64,delimiter='\n')

	evec_scalar = np.nansum(evecs,axis=0)[np.newaxis,:]

	output=open(outfile,"wb")

	(bim, _, G) = read_plink(filename)
	snps = bim['snp']
	del(bim)

	ncols = evecs.shape[0]	

	for counter in range(int(np.ceil(G.shape[0]/chunk_size))):
		if verbose:
			print("Reading {}".format((counter+1)*chunk_size))
		
		labels = snps[counter*chunk_size:(counter+1)*chunk_size]

		genos = G[counter*chunk_size:(counter+1)*chunk_size,:].compute()
		
		p = np.nanmean(genos,axis=1)/2		

		if missing:
			genos = np.nan_to_num(genos)

		scores = np.dot(genos,evecs)

		scores = scores - 2*np.dot(p[:,np.newaxis],evec_scalar)

		scores = scores / np.sqrt(2*p*(1-p))[:,np.newaxis]

		statistic = (1/evals) * (scores**2)

		statistic = np.insert(statistic.astype(str),0,labels,axis=1)
		
		np.savetxt(output,statistic,delimiter="\t",fmt="%s")	

	output.close()
	exit(0)
