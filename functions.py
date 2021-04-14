from Bio import SeqIO, pairwise2, Seq
import operator,re

# theoretical IgV amplicon sequence
q0 = "AAGATCACCTGCTCCGGGGATAGCAGCTACGCTGGAAGTTACTATTATGGCTGGTACCAGCAGAAGGCACCTGGCAGTGCCCCTGTCACTGTGATCTATGACAACACCAACAGACCCTCGAACATCCCTTCACGATTCTCCGGTTCCAAATCCGGCTCCACAGCCACATTAACCATCACTGGGGTCCGAGCCGACGACAATGCTGTCTATTACTGTGCGAGTACAGACAGCAGCAGTACTGCAGTTGGGGCCGGGACAA"

# extracted by hand from Reynaud 1987, Cell
pseudo_changes = ["T13A-C14A-C15T-G16C-G17T-G19A-A20T-G23T-G26A-T28C-G73A-C80T-G91C-T99C-G100T-C102T-A103G-C105T-A106G-C107A-C108T-A109G-C111G-G154A-G163C-C164A-G188A-C195T-A199G-T201G-C218G-G221A-A223G-C224A-A235G-C242G-A243T",
"A1G-A20G-A22T-G23A-T28G-A29G-GGCTGGTACCAGCAGAAGGCA49G-G100A-C102T-A103G-A104G-C107A-C108T-A121G-C147T-A148G-A149G-G163A-C164A-G188A-C195G-A199G-T201G-G205A-C218G-G222C-C224T-G236T-T237C-T240A-G244A",
"C61T-A62G-A64C-A65T-C68T-A69T-A88T-G100A-A101G-A104G-G120A-C153T-G163C-C164A-G188A-C195T-A199G-T201G-C218G-T222C-A235T-G236A-A238G-C239T-C242G-A243T-G244A",
"G18GAGT-T21C-C30TGGT-T48C-G49A-T55C-G67T-A69T-G91C-G100A-A101G-A106G-C107A-C111G-A121G-A148G-A149C-G188A-C191T-C195G-A199G-T201G-G219A-G221A-T222C-A223T-C224T-A225C-A243T",
"A20G-T21TGGTAGC-G67T-A69T-C83T-C102A-A104G-C111G-A149C-C153T-G163A-C164G-G188A-C195G-A199G-T201G-A212T-C218G-T222C-A223T-C224A-A225C-C234T-A235G-G236C-A238G-C239G-G241T-C242A-A243T",
"A4C-T5A-C8T-G16A-G19A-A20T-A22G-G23T-C24T-G26A-T28A-T48C-A56T-G67T-A101G-C102T-A10b4G-A110G-C111T-G120A-A150G-C195G-A199G-T201G-C218G-A223G-C234T-A235G-A238G-C242G-A243T-G244C",
"A1G-G19A-A20G-C24T-A25G-T28A-C29G-C57T-G67T-G100T-A106G-C107A-C111G-A121G-A162G-C164G-G188A-C195G-A199G-T201G-C218G-T222C-C224G-A225G-A235T-G236A-A238G-C239T-C242G",
"A20G-A22G-C30T-A148C-A149T-G163A-C164A-G188A-C191T-C195G-A199G-T201G-G205A-A212T-C218G-T222C-A223T-C224A-A225C-C228A-A229G-G236C-C239G-A243T",
"A1G-A2G-C8G-C9T-T10C-G11A-C12G-G19A-A20G-G26A-C27T-T28G-C30G-A56T-C61G-GGCA66A-C75T-C80G-C81T-C82G-C87G-T97C-T99C-G100A-A101G-C102T-A103G-G134A-G142T-A148C-A149C-C153T-A162G-G163C-C164A-C167A-A168G-T169C-G182A-G188A-C192T-C195G-T201G-C207A-T214C-C218G-G219A-A220G-G221A-T222C-A223T-C224G-G226A-C228G-A229G-G230C-C231A-A232G-G233C-C234A-A235G-G236T-T237G-C239T-T240C-G241T-C242A-A243T-G244T",
"A20G-A22G-C24T-A25G-T28A-A29G-A56T-G67T-A69T-G91C-G100A-C107A-A121G-A150G-C195G-A199G-T201G-C218G-A223G-G236C-A238T-C239A-A243T",
"A2G-G19A-A20G-A22G-C27T-T28A-C30CAGTGCT-A56T-G67T-T72G-C83T-T99C-G100T-A101G-C102G-C105T-A106G-C107A-C111G-G120A-A150G-C195G-A199G-T201G-C218G-A223G-C234T-A235G-A238G-C242G-A243T-G244C",
"A4C-A20G-T21TAGCAGT-A22G-A25T-C27T-T28G-A29G-C30T-C57T-G67T-G91C-G100A-C107T-C111G-A121G-A162G-C164G-G188A-C195G-A199G-T201G-C218G-G221A-A223G-C242A-A243T-G244A",
"A2G-C6A-C12G-C15T-G18A-A20G-A25T-G26A-C27T-T28G-C30T-C51T-G66A-G67T-G91C-G100A-A101G-C107A-C111G-A121G-T145G-A148C-A149T-C195G-A199G-T201G-A212T-C218G-G219T-A220T-G221C-T222C-A223T-C224G-A225G-C228A-C239A-G241C-A243T",
"C15T-A20G-T21C-A25T-G26A-C27T-T28A-A29G-A56T-G67T-A69T-C83T-T99C-G100T-A101G-C102G-A103G-C105T-A106G-C107A-C108T-A109G-C111G-A121G-A122G-A149C-A162G-C164G-C192T-C195G-A199G-T201G-C218G-T222C-A223T-C224T-A225C-G236C-A238G-C239A-G241A-C242G-A243T",
"A1G-A2G-C8G-C9T-T10C-G11A-C15T-G19A-A20G-A25C-G26A-C27T-T28G-C30G-A56T-C61G-GGCA66A-C75T-C80G-C81T-C82G-C87G-T97C-T99C-G100A-A101G-C102T-A103G-G134A-G142A-A149C-A150C",
"C15A-G18C-G19A-A20G-A25G-G26A-C27T-T28G-A29C-C30T-G67T-A69T-A101G-C102T-A109G-A110G-C111T-C116T-C119A-A121G-G134T-T136A-C138G-T139G-C140G-C141G",
"A2G-C8G-C15T-G19A-A20G-C24T-A25G-C27T-T28G-C51CAGTAACTGTCACGGG-T52C-G53A-T55C-G63A-G67A-A69G-C81A-C82T-A101C-C102T-A104G-T118A-G120C-C147T-A148G-A149G-G163A-C164A-C195G-A199G-T201G-C218G-G219T-A220G-T222C-A223T-C22A-A225C-A229G-G236C-A238G-C239A-A243T",
"A1G-A20G-C27CGATGATGGAAGTTAT-G34A-T40C-G67T-A69T-G100A-A101G-C107A-A109C-C111G-A149C-G163A-C164G-G188A-A199G-T201G-G204A-C218G-T222C-A223T-C224G-A225G-C228T-G236C-A238G-C239A-T240A-A243T",
"T13A-C14G-C15T-G18A-A20G-G23T-A25G-T28C-C30GTGG-G67T-C80T-C81G-T97A-A101C-C102T-A132C-A148C-A149T-A150C-G155T-C156A-C158T-G163C-C164A-G188A-C195G-C198T-A199G-T201G-C218G-G219T-A220G-G221C-T222C-A223T-C22A-A225C-A235T-G236A-A238G-C239T-C242G",
"C15T-G18TAGC-G19A-A20G-T21C-T55A-A56T-G66A-G67T-C96T-C102T-A204G-C107G-A110G-A121G-C152T-A162G-G188A-C191T-C195G-A199G-T201G-C218G-G219T-A220G-G221A-T222C-A223T-C224G-A225G-C231G-G236A-C239G-A243T",
"T13A-G16A-A22G-G23A-A29G-C30G-T55C-A59T",
"T13C-C15A-G16C-G17A-G19A-A20C-T21C-G23T-A25T-G26T-C27T-T48C-G49A-T52G-G67A-A69G-C81A-C82T-A101C-C102T-A104G-T118A-G120C-T125C-C127A-T129G-T130G-C131G-A132G",
"A1G-A20G-G26A-T28A-C57T-G67T-C96T-G100C-C102A-C111G-A114G-A121G-A122G-C123T-A124C-C127T-C133T-G134A-T136C-C138T-C141T-G142T-A148T-A149C-C159T-A160G-G163A-C164A-G188A-G193A-C195G-C198T-A199G-T201G-C218G-G219T-A220G-T222C-A223T-C224G-A225G-C228T-A232T-A235T-T237C-A238C-C239T-T240G-G241A-C242G-A243T-G244C",
"A1G-G91C-T144C-A148C-A149C-T151G-T152G-G154A-C156A-T157A-G163A-C164G-G188A-C191T-C195G-A199G-T201G-T216G-C218G-G221A-T222C-A223T-C224G-A225G-C228T-A229G-G230A-G236C-A238G-C239A-A243T",
"A1G-A2T-A4T-T5A-C8B-C15T-G29A-A20G-T21TAAC-G23A-G26A-T28A-G67A-G85C-G100C-C102A-C108A-C111GC117T-G120A-A121G-A122G-C123T-C126A-C127T-C131T-C133T-G134A-A135C-T139A-C140G-G142T-T144G"]

# read in Sanger sequencing results for the pre-expansion clones
starting_seqs = dict()
for record in SeqIO.parse("Starting_clone_Sanger_trimmed_first11bpsremoved.txt", "fasta"):
    starting_seqs[record.id] = record.seq

# change mutation code
def splitChange(string):
    # -1 to get back to the Python-style numbering
    pos = int("".join([str(f) for f in string if f.isdigit()]))-1
    ref = re.split('[0-9]+', string)[0]
    alt = re.split('[0-9]+', string)[1]
    return([ref, pos, alt])

# this script will inject a given mutation set into a sequence
def putChangesIntoSeq(original, changes, verbose = False):
    ret = list(original)
    doit_later = []
    for c in changes:
        ref, pos, alt = splitChange(c)
        if len(ref) == 1 and len(alt) == 1:
            ret[pos] = alt
        else:
            doit_later.append(c)
    push = 0
    for c in doit_later[::-1]:
        ref, pos, alt = splitChange(c)
        del ret[pos:pos+len(ref)]
        ret.insert(pos, list(alt))
    ret = [item for sublist in ret for item in sublist]
    return("".join(ret))    

# reconstruct pseudogene sequences (e.g. inject the psuedogene changes from the paper into the theoretical sequence)
pseudo_sequences = []

pseudo_changes = list(map(lambda x: x.split("-"), pseudo_changes))

for pc in pseudo_changes:
    ret = putChangesIntoSeq(q0, pc)
    pseudo_sequences.append(ret)
    
# generate pseudogene change sets relative to a given pre-expansion clone
def findRelativePseudoChanges(pseudo_sequences, genotype, barcodes):
    ret = []
    x = str(starting_seqs[genotype + "_" + barcodes])
    for ps in pseudo_sequences:
        # pairwise alignment to get the indels right
        pp = pairwise2.align.globalms(x, ps, 2, -0.8, -6, -1)
        consensus = pp[0][0]
        ps_aligned = pp[0][1]
        changes = []
        correction = 0
        for j in range(len(consensus)):
            if ps_aligned[j] != consensus[j]:
                if consensus[j] == "-":
                    correction -= 1
                changes.append(consensus[j] + str(j+1+correction) + ps_aligned[j])
        ret.append(changes)
    # "ret" will contain the pseudogene differences relative to a given genotype's initial IgV sequence
    return(ret)

# produce a sorted contingency table of a sequence set
# optionally only include those that are more abundant than a relative treshold
def countSeqClasses(what, abundance_treshold = 0):
    store = dict()
    for j in what:
        i = str(j)
        try:
            store[str(i)] += 1
        except:
            store[str(i)] = 0
            store[str(i)] += 1
    store_sorted = sorted(store.items(), key = operator.itemgetter(1), reverse = True)
    ret = []
    for i in store_sorted:
        if i[1] > abundance_treshold:
            ret.append(i)
    return(ret)

# find differences in a given amplicon relative to the pre-expansion initial sequence
# mask: region to omit
def findSeqClassChanges(seq, consensus, mask = None):
    pp = pairwise2.align.globalms(seq, consensus, 2, -0.8, -6, -1)
    scl_aligned = pp[0][0]
    cons_aligned = pp[0][1]
    changes = []
    for j in range(len(scl_aligned)):
        if scl_aligned[j] != cons_aligned[j]:
            if mask is not None and j not in mask:
                changes.append(cons_aligned[j] + str(j+1) + scl_aligned[j])
    return(changes)

# categorize the differences in a given amplicon considering the clone-specific pseudogene changes
# GC events: collapse multiple differences into one event
# GC events: only consider as GC if at least one pseudogene contains all diferences
def determineChangeClass2(relativepseudochanges, classchanges):
    rett = []
    for clch in classchanges:
        which_pseudo = [relativepseudochanges.index(i) for i in relativepseudochanges if clch in i]
        if len(which_pseudo) == 0:
            if splitChange(clch)[0] == "-":
                rett.append(["INS", clch])
            elif splitChange(clch)[2] == "-":
                rett.append(["DEL", clch])
            else:
                rett.append(["PM", clch])
        else:
            pos = splitChange(clch)[1]
            if len(classchanges) == 1:
                rett.append(["AMB", clch])
                continue
            if clch == classchanges[0]:
                posn = splitChange(classchanges[classchanges.index(clch) + 1])[1]
                posp = -1000
                which_pseudo_at_next = [relativepseudochanges.index(i) for i in relativepseudochanges if classchanges[classchanges.index(clch) + 1] in i]
                which_pseudo_at_prev = ["x"]
            elif clch == classchanges[-1]:
                posn = 1000
                posp = splitChange(classchanges[classchanges.index(clch) - 1])[1]
                which_pseudo_at_next = ["x"]
                which_pseudo_at_prev = [relativepseudochanges.index(i) for i in relativepseudochanges if classchanges[classchanges.index(clch) - 1] in i] 
            else:
                posn = splitChange(classchanges[classchanges.index(clch) + 1])[1]
                posp = splitChange(classchanges[classchanges.index(clch) - 1])[1]
                which_pseudo_at_next = [relativepseudochanges.index(i) for i in relativepseudochanges if classchanges[classchanges.index(clch) + 1] in i]
                which_pseudo_at_prev = [relativepseudochanges.index(i) for i in relativepseudochanges if classchanges[classchanges.index(clch) - 1] in i]
            if min(abs(pos-posp), abs(pos-posn)) > 9 or len(set(which_pseudo).intersection(which_pseudo_at_next + which_pseudo_at_prev)) == 0:
                rett.append(["AMB", clch])
            else:
                try:
                    if rett[-1][0] == "GC":
                        rett[-1] = ["GC", rett[-1][1] + "_" + clch]
                    else:
                        rett.append(["GC", clch])
                except IndexError:
                    rett.append(["GC", clch]) 
    return(rett)               


# calculate gene conversion tract length
def getGCLength(gc_event, pseudo_changes, stat):
    potential_pseudogenes = []
    y = gc_event.split("_")[1:-1]
    for pch in pseudo_changes:
        if all([a in pch for a in y]):
            potential_pseudogenes.append(pseudo_changes.index(pch))
    if len(potential_pseudogenes) == 0:
        return([-1])
    else:
        ll = []
        for p in potential_pseudogenes:
            s = pseudo_changes[p].index(y[0])
            e = pseudo_changes[p].index(y[-1])
            if e-s > len(y):
                ll.append(-2)
            else:
                if s > 0:
                    ss = pseudo_changes[p][s-1]
                else:
                    ss = 'x0x'
                if e < len(pseudo_changes[p])-1: 
                    ee = pseudo_changes[p][e+1]
                else:
                    ee = "x262x"
                l = int(splitChange(ee)[1]) - int(splitChange(ss)[1])
                #l = int( splitChange(pseudo_changes[p][e])[1]) - int( splitChange(pseudo_changes[p][s])[1])
                ll.append(l)
        
        if all([x < 0 for x in ll]):
            return([-2])
        else:
            ll = [x for x in ll if x > 0]
        
        if stat == "max":
            ll = np.max(ll)
        if stat == "min":
            ll = np.min(ll)
        if stat == "average":
            ll = np.average(ll)
    return([int(splitChange(ss)[1]), ll])

def getBQfromDictname(name):
    bq = name.split("_")[-1]
    bq = int(bq)
    return(bq)
