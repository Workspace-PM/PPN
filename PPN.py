def get_memory_usage():
    process = psutil.Process()
    return process.memory_info().rss

def metric_finding(ss):
	a = []
	data=np.array([2, 3, 5, 7, 2, 3, 7, 5, 2, 5, 3, 7, 2, 5, 7, 3, 2, 7, 3, 5, 2, 7, 5, 3, 3, 2, 5, 7, 3, 2, 7, 5, 3, 5, 2, 7, 3, 5, 7, 2, 3, 7, 2, 5, 3, 7, 5, 2, 5, 2, 3, 7, 5, 2, 7, 3, 5, 3, 2, 7, 5, 3, 7, 2, 5, 7, 2, 3, 5, 7, 3, 2, 7, 2, 3, 5, 7, 2, 5, 3, 7, 3, 2, 5, 7, 3, 5, 2, 7, 5, 2, 3, 7, 5, 3, 2])
	length=len(data)
	k=0
	while(k<96):
		A=data[k]
		C=data[k+1]
		G=data[k+2]
		T=data[k+3]
		k=k+4
		arr=np.array([])
		seq=np.array([])
		n=len(ss)
		letter_to_num = {'A': A, 'C': C, 'G': G, 'T': T}
		# loop over the sequence and replace each letter with its numeric value
		arr = [letter_to_num[letter] for letter in ss]
		#convert the list of numeric values to a string
		#arr = ''.join(str(num) for num in num_seq)
		arr=np.array(arr)
		size=len(arr)
		r=1
		i=0
		t=2  
		l=4
		p=arr[i] *arr[i+1] *arr[i+2] *arr[i+3] *arr[i+4]
		seq = np.append(seq,p)
		i=i+t	
		p=arr[i-2] *arr[i-1] *arr[i] *arr[i+1] *arr[i+2] *arr[i+3] *arr[i+4]
		seq = np.append(seq,p)
		new_array=[]
		for d in range(0,size):
			new_array.append(arr[d])
		# Define window size and step size
		window_size = 9
		step_size = 2
		# Calculate number of subarrays
		num_subarrays = 1 + (len(new_array) - window_size) // step_size

		# Create view of array with window size and step size
		subarrays = np.lib.stride_tricks.as_strided(new_array, shape=(num_subarrays, window_size), strides=(arr.strides[0]*step_size, arr.strides[0]))
		c=num_subarrays
		
		#print(c)
		subarray_products = subarrays.prod(axis=1)
		i=2*c+2
		seq=np.append(seq,subarray_products)	
		i=i+2	
		if(i==(size-4)):
			p= arr[i-4]*arr[i-3] *arr[i-2] *arr[i-1] *arr[i] *arr[i+1]*arr[i+2]*arr[i+3]
			seq = np.append(seq,p)
			i=i+2
		elif(i==(size-3)):
			p= arr[i-4]*arr[i-3] *arr[i-2] *arr[i-1] *arr[i] *arr[i+1]*arr[i+2]
			seq = np.append(seq,p)
			i=i+2
		elif(i==size-1):
			p= arr[i-4]*arr[i-3] *arr[i-2] *arr[i-1] *arr[i]
			seq = np.append(seq,p)
			i=i+2
		elif(i==(size-2)):	
			p= arr[i-4]*arr[i-3] *arr[i-2] *arr[i-1] *arr[i]*a[i+1]
			seq = np.append(seq,p)	
			i=i+2
		value=np.sum(seq)
		a.append(value)
	return(a)

def prime():
	fasta_file = open(sys.argv[1], "r")
	sequences = list(SeqIO.parse(fasta_file, "fasta"))
	# Extract the IDs and sequences into separate lists
	seq_id = [seq.description for seq in sequences]
	seq = [str(seq.seq) for seq in sequences]
	
	print("\n----------Method: PPN----------------")
	print("Total number of species: ",len(seq_id))
	seq_id=np.array(seq_id)
	length=len(seq_id)
	all_vector=[]
	one_d_metric=[]
	for i in range(length):
		s=seq[i]
		s = str([element.upper() for element in s])
		s = re.sub('[^ACGT]', '', s)
		one_d_metric=metric_finding(s)
		all_vector.append(one_d_metric)
		print("Sequence of",seq_id[i],"was processed, Length",len(s))
	vec = len(all_vector)
	
	#Generate all combinations of 2 base addresses
	combinations = list(itertools.combinations(range(vec), 2))
	distance=np.array([])	
	n=len(combinations)
	for i in range(n):
		first_pair = combinations[i]
		element1, element2 = first_pair
		point1=np.array(all_vector[element1])
		point2=np.array(all_vector[element2])
		
		#Euclidean metric
		p = np.linalg.norm(point1 - point2) 
			
		#Manhatten
		#p=np.sum(np.abs(point1-point2))
			
		#cos_sim = np.dot(point1,point2) / (np.linalg.norm(point1) * np.linalg.norm(point2))
		# calculate the cosine distance
		#p = 1 - cos_sim
		distance = np.append(distance,p)
		
	#Distance matrix section
	maximum=max(distance)
	if(maximum!=0):		
		distance=np.around(distance/maximum,3)
	else:
		print("all species are same")
	distance_matrix = squareform(distance)
	filename = sys.argv[2]
	# write the distance matrix to the CSV file
	with open(filename, mode='w', newline='') as file:
		writer = csv.writer(file)
		# write header row
		header_row = [''] + list(seq_id)
		writer.writerow(header_row)
		# write data rows
		for i, row in enumerate(distance_matrix):
			writer.writerow([seq_id[i]] + list(row))
			
	# Tree section
	inputFile=sys.argv[2]
	outputFile=sys.argv[3]
	pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(inputFile), delimiter=",")
	if(sys.argv[4]=="NJ"):
		tree = Nonetree = pdm.nj_tree()
	elif(sys.argv[4]=="UPGMA"):
		tree = Nonetree = pdm.upgma_tree()
	else:
		print("No/Incorrect tree generation algorithm given.")
		sys.exit(1)
	f = open(outputFile, "w")
	f.write(tree.as_string("newick"))
		
if __name__ == '__main__':
	import psutil
	import itertools
	import numpy as np
	import sys
	import re
	import time
	import math	
	import csv
	from scipy.spatial.distance import cosine, pdist, squareform
	from Bio import SeqIO	
	from scipy.cluster import hierarchy
	import pandas as pd
	from numpy import linalg as LA
	import dendropy
	import os
	from scipy.stats import entropy
	import resource
	
	print("Work in progress...........")
	resource.setrlimit(resource.RLIMIT_AS, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))

	#Track the peak memory usage
	peak_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	start_time=time.time()		
	prime()
	end_time=time.time()
	t=end_time -start_time
	print("\n\nTotal Execution time is:", t,"Sec")	
	final_peak_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	
	# Calculate the peak memory consumption
	peak_memory_consumption = final_peak_memory - peak_memory
	peak_memory_consumption=peak_memory_consumption / 1024
	
	# Print the peak memory consumption
	print(f"Peak memory consumption: {peak_memory_consumption} Megabytes")
	with open("time_memory.txt", "w") as file:
	# Write the value to the file
		file.write("Memory:" + str(peak_memory_consumption) + "\t")
		file.write("Time:" + str(t))
		file.write("\n	")
	
	
	
	
	
	
        	
