import matplotlib.pyplot as plt
plt.switch_backend("agg")
from Bio import SeqIO
from Bio import SeqUtils
from collections import Counter
from collections import OrderedDict
import numpy as np

gc_content = []

for record in SeqIO.parse("/Johnny/data/input/Bacteria/E.coli/K12/ucsd_lane_1/ecoli_mda_lane1.fastq",
                          "fastq"):
    num_qualities = record.letter_annotations["phred_quality"]

    assert (len(record.seq) == len(num_qualities))

    aver_quality = sum(num_qualities) / float(len(num_qualities))
    if aver_quality < 35:
        continue

    bad_num_qualities = [x for x in num_qualities if x < 33]
    if len(bad_num_qualities) >= len(num_qualities) * 0.5:
        continue

    gc_count = 0
    high_quality_length = 0

    for i, c in enumerate(record.seq):
        #if num_qualities[i] >= 35:
        if c == 'G' or c == 'C':
		gc_count += 1
	high_quality_length += 1
    gc_content.append(gc_count*100/high_quality_length)

gc_content_countered = OrderedDict(sorted(Counter(gc_content).items(),
                                          key=lambda t: t[0]))

x, y = zip(*gc_content_countered.items())

plt.plot(x, y)
plt.xlabel("%GC")
plt.ylabel("Number of reads")
plt.title("gc content")
plt.savefig("gc_content.png")



