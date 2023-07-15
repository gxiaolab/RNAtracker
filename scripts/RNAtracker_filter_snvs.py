#prepare input files
#annotate the SNVs with gene and context
import os
import pickle
import sys
import json
import numpy as np
import pybedtools
import argparse
from progress.bar import Bar
from scipy.sparse import csr_matrix

def show_value(s):
    """
    Convert unicode to str under Python 2;
    all other values pass through unchanged
    """
    if sys.version_info.major == 2:
        if isinstance(s, unicode):
            return str(s)
    return s

# get the genomic context of SNVs for multiple regions
PRIORITY_REGION_LIST = ["exon", "3UTR", "5UTR", "intron", "noncodingExon", "noncodingIntron"]
PRIORITY_CODES = dict(zip(PRIORITY_REGION_LIST,[1,2,3,44,5,66]))

def get_best(region=None):
    if not isinstance(region, list):
        print("{} is not recognized".format(region))
        exit(0)
    for i in PRIORITY_REGION_LIST:
        if i in region:
            return i
            #if i.startswith("noncoding"):
            #    return "noncoding"
            #else:
            #    return i
    return None

# obtain the genomic context of RDD input, will generate a corresponding BED file and a .annot file
def get_context(ifn=None, annotation=None):
    # chr4 725079 G>A nan 7:5:0 +
    # convert to bed format: $4 name, $5 score, $6 strand
    # output bed file and annotation file:
    ofn1 = opts.prefix + '/' + os.path.basename(ifn) + ".bed"
    ofn2 = opts.prefix + '/' + os.path.basename(ifn) + ".annot"
    if os.path.isfile(ofn2):
        assert os.path.isfile(ofn1)
        return 
    else:
        with open(ofn1, "w") as outputf:
            with open(ifn) as inputf:
                for nline, strline in enumerate(inputf):
                    listline = strline.split()
                    chrom = listline[0]
                    site = listline[1]
                    strand = listline[5]
                    base = listline[2]
                    rcs = listline[4]
                    outputf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, site, site, base, rcs, strand))
        f_snv = pybedtools.BedTool(ofn1)
        f_annot = pybedtools.BedTool(annotation)
        print(ofn1, annotation)
        intersect = f_snv.intersect(f_annot, wo=True, s=True) # -s parameter to require same strandedness
        # extract gene and region for each SNV
        dict_snv_region = {}
        for i, istr in enumerate(intersect):
            # istr: 86556 chr3 69045743 69045743 C>T 9:7:0 - chr3 69044595 69047357 TMF1|ENST00000646304.1|intron:013 0 - 0
            lst_itsc = [show_value(i) for i in istr.fields] 
            coord = lst_itsc[0] + ":" + lst_itsc[1] + ":" + lst_itsc[5]
            mut = lst_itsc[3]
            count = lst_itsc[4]
            context = lst_itsc[9]
            gene = context.split("|")[0]
            region = context.split("|")[2].split(":")[0]
            if region in PRIORITY_REGION_LIST:
                env = gene + "_" + region
                if coord not in dict_snv_region:
                    dict_snv_region[coord] = []
                    dict_snv_region[coord].append(mut)
                    dict_snv_region[coord].append(count)
                    dict_snv_region[coord].append(env)
                else:
                    if env not in dict_snv_region[coord]:
                        dict_snv_region[coord].append(env)
                    else:
                        continue
        with open(ofn2, "w") as of:
            for snv in dict_snv_region:
                chrom = snv.split(":")[0]
                site = snv.split(":")[1]
                strand = snv.split(":")[2]
                mut = dict_snv_region[snv][0]
                count = dict_snv_region[snv][1]
                envs = dict_snv_region[snv][2:]
                if len(envs) == 1:  # check if the SNV has only 1 gene and 1 context
                    gene = envs[0].split("_")[0]
                    reg = envs[0].split("_")[1]
                    # chr4 725079 G>A nan 7:5:0 +
                    if reg in PRIORITY_REGION_LIST:
                        of.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, site, mut, "nan", count, strand, gene, reg))
                else:  # SNV with more than 1 context
                    # print(snv,mut,count,envs)
                    genes = []
                    regs = []
                    for e in envs:
                        genes.append(e.split("_")[0])
                        regs.append(e.split("_")[1])
                    # check if the SNV is only in 1 gene
                    genes_nodup = list(set(genes))
                    if len(genes_nodup) == 1:
                        gene = genes_nodup[0]
                        reg = get_best(regs)
                        if reg in PRIORITY_REGION_LIST:
                            of.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format \
                                     (chrom, site, mut, "nan", count, strand, gene, reg))
                    elif len(genes_nodup) > 1:
                        for g in genes_nodup:
                            subr = []
                            gene = g
                            indices = [i for i, istr in enumerate(genes) if istr == gene]
                            for i in indices:
                                subr.append(regs[i])
                            if len(subr) == 1:
                                reg = subr[0]
                            elif len(subr) > 1:
                                reg = get_best(subr)
                            if reg in PRIORITY_REGION_LIST:
                                of.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, site, mut, "nan", count,
                                                                               strand, gene, reg))

def read_rdd(fn):
    res = {}
    with open(fn) as fp:
        for i, istr in enumerate(fp):
            t = istr.split("\t")
            k = ":".join([t[0], t[1], t[5]])
            t3 = [int(t2) for t2 in t[4].split(":")[:-1]]
            res[k] = {'counts': t3, 'allele': t[2]}
    return res

def read_normalized_rdd(filename_list=None,
                        fn_counts=None,
                        gene_site_matrix=None):
    nreplicates = len(filename_list)
    # get the total number of reads of each sample (and replicates)
    counts_dict = {}
    with open(fn_counts) as fp:
        for i, istr in enumerate(fp):
            t = istr.split()
            counts_dict[t[-2]] = int(t[-1])
    print(counts_dict.keys())
    print(filename_list)
    total_cts = [counts_dict[os.path.basename(t).split('.')[0]] for t in filename_list]
    print(total_cts)
    normalize_factor = [np.mean(total_cts) / total_cts[i] for i in range(nreplicates)]
    raw_data = [read_rdd(fn) for fn in filename_list]

    # get the common SNVs in the k replicates
    raw_keys = [set([k for k in sample]) for sample in raw_data]
    common_keys = raw_keys[0]
    for i in range(1, len(raw_keys)):
        common_keys = common_keys.intersection(raw_keys[i])
    common_keys = list(common_keys)
    common_keys.sort()
    print ("number of common keys: {}".format(len(common_keys)))

    #decide the gene and region of common SNVs in the k replicates (common_keys)
    #load the correspondance between gene and site matrix built
    with open(gene_site_matrix,"rb") as fp:
        t_matrix = pickle.load(fp)
        genes_names = t_matrix['genes_names']
        sites_names = t_matrix['sites_names']
        gene_site_matrix = t_matrix['gene_site_matrix']
    print ("shape of gene_site_matrix: {}".format(gene_site_matrix.shape))

    #to get counts, normalized_counts, indices_in_map, and genomic region (site_region) of common_keys
    indices_in_map = [] #indices of common_keys in gene_site_matrix
    site_region = [] #genomic regions of common_keys to extract from gene_site_matrix
    print(gene_site_matrix.shape, gene_site_matrix.max(),gene_site_matrix.min())
    counts = []
    normalized_counts = []
    dict_site_name_to_int = dict(zip(sites_names,range(len(sites_names))))
    with Bar("Processing all the commons keys",max=len(common_keys)) as bar:
        for k in common_keys:
            t = dict_site_name_to_int[k]
            indices_in_map.append(t)
            ti = gene_site_matrix[t].nonzero()[1] #sparse, use toarray(); rows in gene_site_matrix
            #minimum region code of the site: 1 exon, 2 3UTR, 3 5UTR, 44 intron, 5 noncodingExon, 66 noncodingIntron
            site_region.append(gene_site_matrix[t,ti].min()) #this is not the ultimate region of a site, just used to decide if the site should be kept. eg: 
#chr15	92902553	G>A	nan	9:7:0	+	CHD2	3UTR
#chr15	92902553	G>A	nan	9:7:0	+	AC013394.1	intron 
            #then extract the counts and normalized_counts
            cts = []
            ncts = []
            for n in range(nreplicates):
                c = np.asarray(raw_data[n][k]['counts'])
                s = normalize_factor[n]
                cts.append(c)
                ncts.append(c * s)
            counts.append(cts)
            normalized_counts.append(ncts)
            bar.next()

    selected_gene_site_matrix = gene_site_matrix[indices_in_map]
    counts = np.asarray(counts)  # shape = (nsites,nreplicates,2) #2 for ref and alt
    normalized_counts = np.asarray(normalized_counts)
    sum_ncts = normalized_counts.sum(axis=2)
    allelic_ratio = normalized_counts.max(axis=2) / normalized_counts.sum(axis=2)
    # decide the maximum allowed read coverage for an SNV, to exclude the SNVs with extremely high coverage (and allelic ratio near 1.0 ),
    # which may due to experimental bias.
    maximum_allowed_ncts = np.max(sum_ncts[np.where(np.abs(allelic_ratio - 0.5) < 0.05)]) 
    print("maximum ncts {}".format(np.max(normalized_counts)))
    print("maximum allowed ncts {}".format(maximum_allowed_ncts))
    valid_snv = np.all(sum_ncts < maximum_allowed_ncts, axis=1)
    return {'common': {'keys': common_keys,
                       'indices_in_map': indices_in_map,
                       'counts': counts,
                       'normalized_counts': normalized_counts,
                       'allelic_ratio': allelic_ratio,
                       'valid_snv': valid_snv,
                       'region': site_region, #just used to decide if the site should be kept
                       'selected_gene_site_matrix':selected_gene_site_matrix},
            'gene_names': genes_names,
            'raw': raw_data}


parser = argparse.ArgumentParser(description="Scripts descriptions here")
parser.add_argument('--rdd', dest="rdd", required=True, nargs="+", help='input rdd files')
parser.add_argument('-c', required=False, help='total mapped read number of each sample')
parser.add_argument('-g', required=False, help='gene annotation',
                    default="/u/project/gxxiao3/lingzhan/project/bruseq_encode/annotation/GRCh38_gencode36.bed")
parser.add_argument('--prefix', dest='prefix', help='prefix string for all the outputs')
opts = parser.parse_args()

#first get the genomic context of each input rdd files
for t in opts.rdd:
    dict_snv_context=get_context(ifn=t, annotation=opts.g)

#read the gene-site map to a json file
annotf = opts.prefix + '/' + os.path.basename(opts.rdd[0]) + ".annot" 
print (annotf) #chr4    2837932 G>A     nan     7:6:0   +       SH3BP2  3UTR
sites2genes = {}
raw_site_names = []
raw_gene_names = []
fn_map = "{}/gene_site_map.json".format(opts.prefix)
if not os.path.isfile(fn_map):
    with open(annotf, "r") as fp:
        with Bar("Building the gene sites map") as bar:
            for i, istr in enumerate(fp):
                t = istr.split()
                site = t[0]+":"+t[1]+":"+t[5]
                gene = t[6]
                context = t[7]
                if site in sites2genes:
                    sites2genes[site].append((gene,context))
                else:
                    sites2genes[site] = [(gene,context)]
                raw_site_names.append(site)
                raw_gene_names.append(gene)
                bar.next()
    t = {'sites2genes': sites2genes,
         'raw_site_names':raw_site_names,
         'raw_gene_names':raw_gene_names}
    with open(fn_map, "w") as fp:
        json.dump(t, fp)
else:
    with open(fn_map, "r") as fp:
        t = json.load(fp)
        sites2genes = t['sites2genes']
        raw_gene_names = t['raw_gene_names']
        raw_site_names = t['raw_site_names']
    
print("Gene Site map built")

# get gene-site correspondance matrix and save it to a pickle file
fn_matrix = "{}/gene_site_matrix.pickle".format(opts.prefix)
if os.path.isfile(fn_matrix):
    with open(fn_matrix,"rb") as fp:
        t_matrix = pickle.load(fp)
        genes_names = t_matrix['genes_names']
        sites_names = t_matrix['sites_names']
        gene_site_matrix = csr_matrix.todense(t_matrix['gene_site_matrix'])
else:
    genes_names = list(set(raw_gene_names))
    genes_names.sort()
    dict_g = dict(zip(genes_names,range(len(genes_names))))
    sites_names= list(set(raw_site_names))
    sites_names.sort()
    dict_s = dict(zip(sites_names,range(len(sites_names))))
    print("dictionaries are built")
    gene_site_matrix = np.zeros(shape=(len(sites_names),len(genes_names)),dtype=np.int8)
    for sn,values in sites2genes.items():
        for gn,cx in values:
            gene_site_matrix[dict_s[sn],dict_g[gn]] = PRIORITY_CODES.get(cx,0)

    t_matrix = {'genes_names':genes_names,
                'sites_names':sites_names,
                'gene_site_matrix':csr_matrix(gene_site_matrix, dtype=np.int8)}
    with open(fn_matrix,"wb") as fp:
        pickle.dump(t_matrix,fp)

#start to read rdd data of the k replicates
list_rdd_annot=[]
for rdd in opts.rdd:
    rdd_item = opts.prefix + '/' + os.path.basename(rdd) + ".annot" 
    list_rdd_annot.append(rdd_item)

data_all = read_normalized_rdd(filename_list=list_rdd_annot, fn_counts=opts.c, gene_site_matrix=fn_matrix)
fn_alldata = "{}/data_all.pickle".format(opts.prefix)
with open(fn_alldata, "wb") as fp:
    pickle.dump(data_all, fp)
print("All the data are read, and written in {} for plotting and analysis".format(fn_alldata))

