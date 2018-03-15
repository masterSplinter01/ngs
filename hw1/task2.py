from __future__ import division
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from Bio import SeqIO
import itertools

strs = []
count = 0

for record in SeqIO.parse("/Johnny/data/input/Bacteria/E.coli/K12/ucsd_lane_1/ecoli_mda_lane1.fastq",  "fastq"):
    	strs.append(record.seq)
	count+=1
pos = itertools.izip_longest(*strs)
probs = []

for item in pos:
    error_count = 0
    length = 0
    for j in item:
        if j is not None:
            if j == 'N':
                error_count += 1
            length += 1
    prob = error_count/length
    probs.append(prob)

indexes = range(0, len(probs))

plt.plot(indexes, probs)
plt.xlabel("index")
plt.ylabel("probability")
plt.savefig("prob.png")
