import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO, pairwise2
import operator
import os, glob
import pickle
import pyfastx
from functions import *


# each fastq pair belongs to a genotype
# internal barcodes refer to parallel clones
# preprocessing: trimming with Trimmomatic, merging with FLASH2
files=sorted(glob.glob("merged/*extendedFrags.fastq.gz"))
codes=list(map(lambda x: "IgV_" + str(x + 1), [0, 1, 2, 3, 4, 6]))
genotypes=['WT', 'BRCA1', '53BP1', 'BRCA1_53BP1', 'KU70', 'REV1']

# internal barcode-related information
forw_bc = ["TGACCA", "TTAGGC", "ATCACG", "ACGGTT", "CGATGT", "CAGTAC", "GCTATG", "CTCAGT", "AGCGTA", "GAGTCA"]
rev_bc = ["TGGTCA", "GCCTAA", "CGTGAT", "AACCGT", "ACATCG", "GTACTG", "CATAGC", "ACTGAG", "TACGCT", "TGACTC"]
forw_names = ["F" + str(i) for i in range(2, 12)]
rev_names = ["R" + str(i) for i in range(2, 12)]
good_ones = ["F"+str(i)+"R"+str(i) for i in [3,  7, 9]]


# separate each file according to the barcodes
# write amplicons with valid combinations to disk
ee = []
oo = []
os.makedirs("obj/")
for j in [0, 1, 2, 3, 4, 6]:
    fq = pyfastx.Fastq('merged/IgV_' + str(j+1) + '.extendedFrags.fastq.gz')

    dd = dict()
    cc = dict()

    for g in good_ones:
        dd[g] = []
    
    for a in forw_names:
        for b in rev_names:
            cc[a+b] = 0


    end_is_not_valid_barcode = 0
    other_barcode_group = 0

    for i in range(len(fq)):
        s = fq[i].seq
        f=s[2:8]
        r=s[len(s)-8:len(s)-2]
        try:
            address = forw_names[forw_bc.index(f)] + rev_names[rev_bc.index(r)]
            dd[address].append(fq[i].id)
            cc[address] += 1
        except ValueError:
            end_is_not_valid_barcode += 1
        except KeyError:
            other_barcode_group += 1
        if i % 10000 == 0:
            print(str(i) + ":" + " (Not barcode) " + str(end_is_not_valid_barcode) \
                  + " (Bad combination) " + str(other_barcode_group) \
                  + " (Good combination) " + str(sum([len(a) for a in list(dd.values())])))

    ee.append(end_is_not_valid_barcode)
    oo.append(other_barcode_group)

    with open("obj/IgV_" + str(j+1) + "_pyfastqx_all_good_indeces.pkl", "wb+") as f:
        pickle.dump(dd, f, pickle.HIGHEST_PROTOCOL)

        
# find differences between each amplicon and the respective pre-expansion sequence for the given clone
# enumerate the differences, and save to disk for each genotype and clone

for code in codes:
    for good_one in good_ones:
        print("====== " + code + "_" + good_one + " ======") 
        fq = pyfastx.Fastq('merged_notrim/' + code +  '.extendedFrags.fastq.gz')

        with open("obj/" + code + "_pyfastqx_all_good_indeces.pkl", "rb") as f:
            dd = pickle.load(f)

        relative_pseudo_changes = findRelativePseudoChanges(pseudo_sequences, genotypes[codes.index(code)], good_one)

        ddq = {"Unchanged": 0}

        not_included = []
        has_indel = []

        for i in dd[good_one]:
            # class Fastq's id property starts at 1, the python ordering starts at 0
            i -= 1
            if len(fq[i]) == len(starting_seqs[genotypes[codes.index(code)] + "_" + good_one]) + 27:
                if fq[i].seq.find("GAGAAACCGTC") == 8:
                    s = fq[i].seq[19:-8]
                    q = fq[i].quali[19:-8]
                elif fq[i].seq.find("TTGTCCCGGC") == 8:
                    s = Seq.Seq(fq[i].seq)
                    s = str(s.reverse_complement()[19:-8])
                    q = fq[i].quali[::-1][19:-8]
                else:
                    not_included.append(i)
                    continue
            else:
                has_indel.append(i)
                continue
            class_changes = findSeqClassChanges(s, str(starting_seqs[genotypes[codes.index(code)] + "_" + good_one]), mask = [100000])
            change_type = determineChangeClass2(relative_pseudo_changes, class_changes)
            if len(change_type) > 0:
                for y in change_type:
                    if "-" not in y[1]:
                        y = y + [str(q[splitChange(y[1].split("_")[0])[1]])]
                    try: 
                        ddq["_".join(y)] += 1
                    except KeyError:
                        ddq["_".join(y)] = 1
            else:
                ddq["Unchanged"] += 1
            if dd[good_one].index(i+1) % 10000 == 0:
                print(str(dd[good_one].index(i+1)) + " / " + str(len(dd[good_one])) + " / " + str(sum(ddq.values())))
        print("=======================================")        
        # save the ddq object
        with open("obj/" + code + "_" + good_one + "_pyfastqx_all_event_counts.pkl", "wb+") as f:
            pickle.dump(ddq, f, pickle.HIGHEST_PROTOCOL)
            
# categorize into GC, PM and AMB
# base qualities are also taken into account: at least 30 is needed 
# point mutations are split by mutation type
outcomes = open("outcomes_all_noindels.txt", "a")
pm_types = open("pm_types_all_noindels.txt", "a")
stats = open("stats_all_noindels.txt", "a")

outcomes.write("\t".join(["Genotype", "Barcode", "AMB", "GC", "PM"]) + "\n")
pm_types.write("\t".join(["Genotype", "Barcode"] + [x + ">" + y for x in bb for y in bb if x != y]) + "\n")
stats.write("\t".join(["Genotype", "Barcode", "Tested_seqs", "Invalid_barcode", "Invalid_combination", "Good_barcode",
                      "All_events", "Unchanged", "AMB_low_BQ", 'GC_low_BQ', "PM_low_BQ"]) + "\n")

for code in codes:
    
    with open("obj/" + code + "_pyfastqx_all_good_indeces.pkl", "rb") as f:
            dd = pickle.load(f)
    fq = pyfastx.Fastq('merged_notrim/' + code +  '.extendedFrags.fastq.gz')
    
    for good_one in good_ones:
        with open("obj/" + code + "_" + good_one + "_pyfastqx_all_event_counts.pkl", "rb") as f:
            ddq = pickle.load(f)
        df = pd.DataFrame(ddq.items(), columns = ["ID", "Count"])
        new = df["ID"].str.split("_", n = 1, expand = True) 
        df["Category"] = new[0]
        df["ID2"] = new[1]
        vv = [d for d in df.ID2.str.split("_")]
        uu = []
        for v in vv:
            if v is None:
                uu.append(v)
            elif len(v) > 1:
                uu.append(int(v[-1]))
            else: 
                uu.append(-1)
        df["BQ"] = uu

        df["Event"] = df["ID2"].str[:-3]
        df.drop(columns = ["ID", "ID2"], inplace = True)
        df_pm = df[df["Category"] == "PM"]
        df_pm["Muttype"] = df_pm["Event"].str[0] + ">" + df_pm["Event"].str[-1]
        outcomes.write("\t".join([genotypes[codes.index(code)], good_ones[good_ones.index(good_one)]] + \
                        [str(x[0]) for x in df[df["BQ"] >= 30].groupby('Category').sum().values.tolist()]) + "\n")
        pm_types.write("\t".join([genotypes[codes.index(code)], good_ones[good_ones.index(good_one)]] + \
                        [str(x[0]) for x in df_pm[df_pm["BQ"] >= 30].groupby('Muttype').sum().values.tolist()]) + "\n")
        stats.write("\t".join([genotypes[codes.index(code)], good_ones[good_ones.index(good_one)], 
                               str(len(fq)), str(ee[codes.index(code)]) , str(oo[codes.index(code)]), str(len(dd[good_one])), 
                               str(df.Count.sum()), str(ddq["Unchanged"]),
                              str(df[df.BQ < 30][df.Category == "AMB"].Count.sum()),
                              str(df[df.BQ < 30][df.Category == "GC"].Count.sum()),
                              str(df[df.BQ < 30][df.Category == "PM"].Count.sum()) + "\n"]))
        
        print(code + "_" + good_one)

outcomes.close()
pm_types.close()
stats.close()