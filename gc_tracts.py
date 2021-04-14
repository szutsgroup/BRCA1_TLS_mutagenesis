import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyfastx
import pickle
from functions import *
from scipy.stats import  mannwhitneyu

codes=list(map(lambda x: "IgV_" + str(x + 1), [0,1,2,3,4,6]))
genotypes=['WT', 'BRCA1', '53BP1', 'BRCA1_53BP1', 'KU70', 'REV1']
forw_bc = ["TGACCA", "TTAGGC", "ATCACG", "ACGGTT", "CGATGT", "CAGTAC", "GCTATG", "CTCAGT", "AGCGTA", "GAGTCA"]
rev_bc = ["TGGTCA", "GCCTAA", "CGTGAT", "AACCGT", "ACATCG", "GTACTG", "CATAGC", "ACTGAG", "TACGCT", "TGACTC"]
forw_names = ["F" + str(i) for i in range(2, 12)]
rev_names = ["R" + str(i) for i in range(2, 12)]
good_ones = ["F"+str(i)+"R"+str(i) for i in [3, 7, 9]]



# calculate and visualize GCtract lengths
plt.figure(figsize = (10, 10))

colors = ["red", "blue", "green", "orange"]
code_subset = ["IgV_1", "IgV_2", "IgV_3", "IgV_4"]
xpos = 0
ll = []
ss = []

for code in code_subset:
    oo = []
    xpos += 1
    ll.append(genotypes[codes.index(code)])
    
    for g in good_ones:
        with open("obj/" + code + "_" + g + "_pyfastqx_all_event_counts.pkl", "rb") as f:
            ddq = pickle.load(f)
        gcs = {k: ddq[k] for k in ddq if "GC" in k if getBQfromDictname(k) >= 30}
        pseudo_changes = findRelativePseudoChanges(pseudo_sequences, genotypes[codes.index(code)], g)
        for k in gcs.keys():
            oo.append(getGCLength(k, pseudo_changes, "min") + [gcs[k]])


    x = [a[1] for a in oo if len(a) == 3]
    y = [a[2] for a in oo if len(a) == 3]
    xx = np.repeat(x, y)
    ss.append(xx)
    
    plt.scatter(x = zz, y = xx, c = colors[code_subset.index(code)])
    plt.plot([xpos - .5, xpos + .5], [np.median(xx), np.median(xx)], linewidth = 6, color= "black")
    plt.plot([xpos - .7, xpos + .7], [np.percentile(xx, 25), np.percentile(xx, 25)], linewidth = 4, color= "black")
    plt.plot([xpos - .7, xpos + .7], [np.percentile(xx, 75), np.percentile(xx, 75)], linewidth = 4, color= "black")
    plt.text(xpos + .6, np.median(xx), round(np.median(xx)), fontsize = 12)

    xpos += 1        

plt.xticks([1, 3, 5, 7], ll)
top = max(sum([i.tolist() for i in ss], []))
sss = [np.bincount(a) for a in ss]
s4 = []
for a in sss:
    tmp = np.zeros(top+1)
    tmp[:len(a)] = a
    s4.append(tmp.tolist())

for i in [1, 2, 3]:
    mwu = mannwhitneyu(x = s4[0], y =s4[i])
    x1, x2 = 1, 1+(i*2)  
    y, h, col = 130 + (i-1)*5, 1, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, str(round(mwu[1], 8)), ha='center', va='bottom', color=col)

plt.ylabel(" Frequency weighted GC tract length")
plt.show()
