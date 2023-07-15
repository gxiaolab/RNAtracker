#############################################################
####### To identify allele-specific expression of gene ######
####### by Ling Zhang                                  ######
####### 2021-04-21                                     ######
#############################################################

# part1: identify ASE SNVs using similar strategy to BEAPR
# part2: identify ASE genes based on all the SNVs in that gene

import os
import sys
import json
import numpy as np
from scipy.sparse import csr_matrix
import scipy 
from statsmodels.stats.multitest import multipletests
import pybedtools
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from scipy import special
import matplotlib.pyplot as plt
import argparse
import joblib
import pickle
from progress.bar import Bar
from libloess_ar import Myloess, plot_loess,test_new_site, test_new_sites
from filtergenes import filtergenes,filtergenes_sparse2
import logging
from datetime import datetime

parser = argparse.ArgumentParser(description="Scripts descriptions here")
parser.add_argument('--in-prefix', dest='in_prefix', help='prefix string for all the inputs') 
parser.add_argument('--prefix', dest='prefix', help='prefix string for all the outputs')
parser.add_argument('--exclude-intronic', dest="ei", action='store_true', help='if only use SNVs in non-intronic regions')
parser.add_argument('--filter-test-set',dest='ftest',nargs=2,type=int, default = [-1,-1],help='[n1,n2],n1 is the minimum value for reads, n2 is the minimal values for total reads')
parser.add_argument('--filter-regression-set',dest='freg',nargs=2,type=int, default = [-1,-1],help='[n1,n2],n1 is the minimum value for reads, n2 is the minimal values for total reads')
parser.add_argument('--debug',dest='debug',action='store_true',help='debug mode, more outputs')
parser.add_argument('--plot-loess',dest='plot',action='store_true',help='plot the results of loess')
parser.add_argument('--plot-sigma',dest='plotsigma',action='store_true',help='plot the distribution of the sigma of delta allelir ratio in control and testable genes')
parser.add_argument('--regression',dest='r',default="cv2vslog2mu",help='regression x and y',choices=["cv2vslog2mu","cvvslog2mu"])
parser.add_argument('--graph-labels',dest="glabels", nargs=2, default=["log2mu",'cv2'], help="graph labels for regression, x and y")
parser.add_argument('--use-allelic-ratio',dest="usear", action='store_true', help="if use the sigma of allelic ratios as y in regression")
parser.add_argument('--adjustP',dest="ap", default="None", help="FDR correction function in multipletests, fdr_bh,fdr_tsbh,fdr_by,hs, if None, no adjust used")
parser.add_argument('-p', default=0.05, type=float, help='cutoff for Pvalue to determine ASE SNVs')
parser.add_argument('--nsigma',dest="nsigma", default=1.645, type=float, help="SNVs within allelic ratio within nsigma is considered as ASE") 
opts = parser.parse_args()

outall_prefix = "{}/T{}_M{}_ei{}_AR{}_nsigma{}_p{}_ap{}".format(opts.prefix,opts.ftest[1], opts.ftest[0],opts.ei,opts.usear,opts.nsigma,opts.p,opts.ap)

logger = logging.getLogger('ASE')
now = datetime.now()
fh = logging.FileHandler('{}_ASE_run.log'.format(outall_prefix))
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)
if opts.debug:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

logger.critical("INPUT options = {}".format(opts))
logger.debug("the input parameters are {}".format(opts))

#read common SNV of the k replicates
fn_dataall = "{}/data_all.pickle".format(opts.in_prefix) 
with open(fn_dataall, "rb") as fp:
    data_all = pickle.load(fp)
logger.info("All the data are loaded from {}".format(fn_dataall))

logger.info("Start to select SNVs for regression and testable SNVs")
#retrieve the part for common SNVs in the k replicates
data_common = data_all['common']
rcounts = data_common['counts']
logger.critical("Number of common SNVs from all replicates: {}".format(len(data_common['keys'])))

c1 = data_common['valid_snv'] #exclude these with extremely high coverage 
#SNVs for regression: should include all valid SNVs (no matter read coverages) when doing SNV level regression
cds_testable = [c1]
cds_regression = [c1]
# for regression
c2rm = rcounts.min(axis=1).min(axis=1) >= opts.freg[0]
c2rt = rcounts.sum(axis=2).min(axis=1) >= opts.freg[1]
c2r = np.all([c2rm,c2rt],axis=0)

c2tm = rcounts.min(axis=1).min(axis=1) >= opts.ftest[0]
c2tt = rcounts.sum(axis=2).min(axis=1) >= opts.ftest[1]
c2t = np.all([c2tm,c2tt],axis=0)

cds_testable.append(c2t)
cds_regression.append(c2r)

if opts.ei:
    c4 = np.asarray(data_common['region']) <10
    logger.critical("Exclude SNVs in intronic and noncodingIntron regions")
    cds_regression.append(c4)
    cds_testable.append(c4)

logger.critical("Exclude SNVs with with allelic ratio >0.95")
allelic_ratios = rcounts[:,:,0]/rcounts.sum(axis=2)

AR_filter1 = allelic_ratios.max(axis=1)<=0.95
AR_filter2 = allelic_ratios.min(axis=1)>=0.05
true_counts = 0
false_counts = 0
for item in AR_filter1:
    if item == True:
        true_counts += 1
    elif item == False:
        false_counts+= 1

logger.critical("Number of SNPs filtered out with AR > 0.95: " + str(false_counts))

true_counts = 0
false_counts = 0
for item in AR_filter2:
    if item == True:
        true_counts += 1
    elif item == False:
        false_counts+= 1

logger.critical("Number of SNPs filtered out with AR < 0.05: " + str(false_counts))
c5 = AR_filter1 & AR_filter2

cds_regression.append(c5)
#better as a choice
cds_testable.append(c5)

selecting_regression = np.where(np.all(cds_regression,axis=0))[0]
logger.critical("Number of sites for regression is {}".format(len(selecting_regression)))
selecting_testable = np.where(np.all(cds_testable,axis=0))[0]
logger.critical("Number of sites for test is {}".format(len(selecting_testable)))

logger.critical("Exclude invalid SNVs (unreasonable read coverages)")
logger.critical("Create SNV sets for regression and test")
data_regression = {}
data_testable = {}
for k, v in data_common.items():
    if k == 'valid_snv':
        continue
    if isinstance(v,list):
        v = np.asarray(v)
    data_regression[k]=v[selecting_regression]
    data_testable[k]=v[selecting_testable]
    if k == 'normalized_counts':
        logger.info("maximum {} is {}".format(k,np.max(data_regression[k])))

data_all['common_valid'] = data_regression
logger.critical("Number of valid SNVs for regression: {}".format(len(data_regression['keys'])))

data_all['common_testable'] = data_testable
logger.critical("Number of testable SNVs under cutoff T={} or M ={}: {}".format(opts.ftest[1], opts.ftest[0],len(data_testable['keys'])))

logger.info("Calculate mu and sigma of major allele counts and allelic ratios of all valid SNVs")
t = data_regression['normalized_counts']
if opts.usear:
    #mu and sigma of allelic ratios of the replicates
    allelic_ratios = np.zeros_like(t) 
    allelic_ratios[:,:,0] = t[:,:,0]/t.sum(axis=2)
    allelic_ratios[:,:,1] = t[:,:,1]/t.sum(axis=2)

    # reshape the allelic_ratios as (nsites, (nreps x 2 ))
    allelic_ratios=allelic_ratios.reshape((t.shape[0],(t.shape[1]*t.shape[2])))
    ratio_sigmas = np.std(allelic_ratios,axis=1)
    sample_mus = t.max(axis=2).mean(axis=1)
    sample_sigmas = ratio_sigmas
    logger.info("Start regression with mu (read coverage of major alleles) and sigma (of allelic ratios) of valid SNVs")
else:
    #mu and sigma of read coverage of major alleles
    coverage_mus = t.max(axis=2).mean(axis=1) #length: nsites, average read coverage of major alleles for each site
    coverage_sigmas =  t.max(axis=2).std(axis=1)
    sample_mus = coverage_mus
    sample_sigmas = coverage_sigmas
    logger.info("Start regression with mu and sigma of read coverage of major alleles of valid SNVs")
loess = Myloess(equation=opts.r, prefix=opts.prefix)
loess.fit_loess(sample_mus,sample_sigmas)

if opts.plot:
   plot_loess(loess, prefix=outall_prefix,labels=opts.glabels)

logger.info("Start test each testable SNVs")
#average coverage of major and minor alleles
ncts = data_testable['normalized_counts']
XMs = ncts.max(axis=2)
Xms = ncts.min(axis=2)
data_testable['XM']=np.mean(XMs,axis=1)
data_testable['Xm']=np.mean(Xms,axis=1)

if opts.usear:
     # find out the sigmas of allielc_ratio
     sigmas = loess.findy(XMs.mean(axis=1)) 
     allelic_ratios = np.zeros_like(ncts) 
     allelic_ratios[:,:,0] = ncts[:,:,0]/ncts.sum(axis=2)
     allelic_ratios[:,:,1] = ncts[:,:,1]/ncts.sum(axis=2)
     #allelic ratios of minor reads
     avg_minor_allelic_ratios = allelic_ratios.min(axis=2).mean(axis=1)
     #allelic ratios of major reads
     avg_major_allelic_ratios = allelic_ratios.max(axis=2).mean(axis=1)
     pvalues = special.erf((avg_minor_allelic_ratios-0.5)/sigmas/np.sqrt(2))+1.0
     data_testable['expected_sigma']=sigmas

else:
    sigmas,pvalues = test_new_sites(loess,XMs,Xms)
    data_testable['expected_sigma']=sigmas

#mutlitple test correction using multiple methods #hs=holm-sidak
methods = ['fdr_bh', 'fdr_tsbh', 'fdr_by', 'fdr_tsbky','fdr_gbs','holm-sidak', 'holm', 'simes-hochberg', 'hommel']

adjust_P_values_set = {}
for met in methods:
    padjs = multipletests(np.asarray(pvalues), method=met, is_sorted=False, returnsorted=False)[1]
    adjust_P_values_set[met]=padjs
adjust_P_values_set['origin']=pvalues
with open("{}_padj.pickle".format(outall_prefix),"wb") as fp:
    pickle.dump(adjust_P_values_set,fp)

logger.debug("number of nan in pvalues {}".format(np.sum(np.isnan(pvalues))))
logger.debug("size of pvalues {}".format(len(pvalues)))
logger.debug("pvalues.max = {}, pvalues.min={}".format(pvalues.max(),pvalues.min()))
logger.critical('number of pvalues smaller than {:.3f} is {}'.format(opts.p, len(np.where(pvalues<opts.p)[0])))
data_testable['original_P'] = np.asarray(pvalues)
if opts.ap == "None":
    data_testable['adjusted_P']=pvalues 
else:
    logger.info("Calculating adjusted p values using {}".format(opts.ap))
    padjs = multipletests(np.asarray(pvalues), method=opts.ap, is_sorted=False, returnsorted=False)[1]
    logger.debug("number of nan in padjs {}".format(np.sum(np.isnan(padjs))))
    data_testable['adjusted_P'] = padjs
    logger.critical('number of adjusted pvalues smaller than {:.3f} is {}'.format(opts.p, len(np.where(padjs<opts.p)[0])))

fn_csdata = "{}/data_selected.pickle".format(opts.prefix)
with open(fn_csdata, "wb") as fp:
    pickle.dump(data_all, fp)
logger.info("Store the data to {}".format(fn_csdata))

logger.info("Start to select control and testable genes:")
gene_site_matrix = data_testable['selected_gene_site_matrix']
if opts.ei:
    #to exclude the intronic and noncodingIntron sites in a gene, mask by coverting 44 and 66 to 0, not consider them in the downstream analysis
    #special cases: a site is in more than 1 genes, in gene1 it is non-intronic, but in gene2 it is intronic, need to exclude the intronic SNVs in gene2
    #example: exclude site chr15:76215497:- in gene TMEM266
    #lingzhan@login3:~/nobackup-gxxiao3/project/bruseq_encode/HepG2/asarp/six_hr/ase$ grep 76215497 ../*.annot
    #../ENCLB378PJW.rdd.final.all_snvs.ENCLB378PJW.rdd.annot:chr15	76215497	C>G	nan	12:8:0	-	ETFA	3UTR
    #../ENCLB378PJW.rdd.final.all_snvs.ENCLB378PJW.rdd.annot:chr15	76215497	C>G	nan	12:8:0	-	TMEM266	intron
    #../ENCLB986YLA.rdd.final.all_snvs.ENCLB986YLA.rdd.annot:chr15	76215497	C>G	nan	18:7:0	-	ETFA	3UTR
    #../ENCLB986YLA.rdd.final.all_snvs.ENCLB986YLA.rdd.annot:chr15	76215497	C>G	nan	18:7:0	-	TMEM266	intron
    exclude = gene_site_matrix>10 
    exclude_index = scipy.sparse.find(exclude)
    gene_site_matrix[exclude_index[0],exclude_index[1]]=0

genes_names = data_all['gene_names'] 
#exclude genes with no site
t = gene_site_matrix.toarray().sum(axis=0)
t2 = np.where(t>0)[0]
gene_site_matrix = gene_site_matrix[:,t2]
genes_names = np.asarray(genes_names)[t2]
logger.debug("number of common keys {}".format(len(data_testable['keys'])))

#cutoff to determine significant SNVs (ASE SNVs)
P_cutoff = opts.p
res_filtered_genes = filtergenes_sparse2(gene_site_matrix, genes_names, data_testable['adjusted_P'], P_cutoff=P_cutoff)

control_genes = genes_names[res_filtered_genes['index_control_genes']]
control_matrix = gene_site_matrix[:,res_filtered_genes['index_control_genes']]
test_genes = genes_names[res_filtered_genes['index_test_genes']]
test_matrix = gene_site_matrix[:,res_filtered_genes['index_test_genes']]
ase_genes = genes_names[res_filtered_genes['index_ase_genes']]
ase_matrix = gene_site_matrix[:,res_filtered_genes['index_ase_genes']]
partase_genes = genes_names[res_filtered_genes['index_part_ase_genes']]
partase_matrix = gene_site_matrix[:,res_filtered_genes['index_part_ase_genes']]
powerful_genes = genes_names[res_filtered_genes['index_powerful_genes']]
powerful_matrix= gene_site_matrix[:,res_filtered_genes['index_powerful_genes']]

logger.info("Control and test genes are selected.")
logger.critical("Number of control genes: {}".format(len(control_genes)))
logger.critical("Number of testable genes: {}".format(len(test_genes)))
logger.critical('number of ASE genes {} (old way)'.format(len(ase_genes)))
logger.critical('number of ASE-part genes {}'.format(len(partase_genes)))
logger.critical('number of powerful genes {}'.format(len(powerful_genes)))


plotting_data = {'control_set_mean':[],
                 'control_set_sigma':[],
                 'test_set_mean':[],
                 'test_set_sigma':[],
                 'ase_set_mean':[],
                 'ase_set_sigma':[],
                 'partase_set_mean':[],
                 'partase_set_sigma':[]}

#with open("{}_ASE_original".format(outall_prefix), "w") as of:
#    for gene in ase_genes:
#        of.write("{}\n".format(gene))

gene_level_mean = []
gene_level_sigma = []

#regression X: average read coverage of all testable SNVs (all replicates)
#regression Y: sigma of delta allelic ratio 
arrayed_control_matrix=control_matrix.toarray()

with Bar("I am working on Loess with control genes, wait",max=arrayed_control_matrix.shape[1]) as bar:
    for ig in range(arrayed_control_matrix.shape[1]):
        t_index = np.where(arrayed_control_matrix[:,ig] > 0 )[0]
        t_norm_cnts = data_testable['normalized_counts'][t_index]
        summed_norm_cnts = t_norm_cnts.sum(axis=1) #shape nsites (in the gene) X 2 , sum the replicates
        delta = summed_norm_cnts.max(axis=1) / summed_norm_cnts.sum(axis=1) - 0.5
        plotting_data['control_set_mean'].append(delta.mean())
        t0 = np.asarray([delta, -1.0 * delta])
        t = t0.flatten()
        sigma = np.std(t)
        plotting_data['control_set_sigma'].append(sigma)
        gene_level_mean.append(summed_norm_cnts.mean())
        gene_level_sigma.append(sigma)
        bar.next()
logger.info("Gene background computed")

gene_level_sigma = np.asarray(gene_level_sigma)
gene_level_mean = np.asarray(gene_level_mean)
logger.debug("gene_level_sigma max = {} and min = {}".format(gene_level_sigma.max(), gene_level_sigma.min()))
logger.debug("gene_level_mean max = {} and min = {}".format(gene_level_mean.max(), gene_level_mean.min()))

loess_gene = Myloess(equation=opts.r, use_interpolation=True, prefix=opts.prefix)
loess_gene.fit_loess(gene_level_mean, gene_level_sigma)
if opts.plot:
    plot_loess(loess_gene, prefix="{}_gene".format(outall_prefix),labels=opts.glabels)

PRIORITY_REGION_LIST = ["exon", "3UTR", "5UTR", "intron", "noncodingExon", "noncodingIntron"]
PRIORITY_CODES = dict(zip([1,2,3,44,5,66],PRIORITY_REGION_LIST))

logger.info ("Start testing each testable gene")
arrayed_test_matrix = test_matrix.toarray()
fn_output = "{}_ASE".format(outall_prefix)
foutput = open(fn_output,"w")
N_ASE = 0
gene_level_sigma_testable = []
gene_level_mean_testable = []
foutput.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format ("Gene", "ASEornonASE", "#SNV", "#SNV-SIG", "#SNV-AR", "Coverage", "Expected_Sigma", "Allelic_Ratio (by ASE SNVs)", "Allelic_Ratio (by testable SNVs)","SNVs"))
with Bar("I am working on Loess with testable genes, wait",max=arrayed_test_matrix.shape[1]) as bar:
    for ig in range(arrayed_test_matrix.shape[1]):
        t_index = np.where(arrayed_test_matrix[:,ig] > 0 )[0]
        t_norm_cnts = data_testable['normalized_counts'][t_index]
        t_raw_cnts = data_testable['counts'][t_index]
        ##way1 for average allelic ratio of the gene: using only the ASE SNVs
        ##extract adjusted P values of each SNV in the gene and then extract the sites with significant adjusted P values
        psig_index = np.where(data_testable["adjusted_P"][t_index]<opts.p)[0]
        sig_norm_cnts =t_norm_cnts[psig_index]
        summed_sig_norm_cnts = sig_norm_cnts.sum(axis=1)
        ars_sig = summed_sig_norm_cnts.max(axis=1) / summed_sig_norm_cnts.sum(axis=1)
        avg_ar_sig = np.mean(ars_sig)
        #way2 for average allelic ratio of the gene: using only all the testable SNVs
        #calculate delta allelic ratio of all testable SNVs in the gene
        summed_norm_cnts = t_norm_cnts.sum(axis=1) #shape nsites (in the gene) X 2 , sum the replicates
        REF_ars = summed_norm_cnts[:,0] / summed_norm_cnts.sum(axis=1)
        ars = summed_norm_cnts.max(axis=1) / summed_norm_cnts.sum(axis=1)
        avg_ar = np.mean(ars)
        delta = ars - avg_ar
        plotting_data['test_set_mean'].append(delta.mean())
        t = np.asarray([delta, -1.0*delta])
        plotting_data['test_set_sigma'].append(np.std(t.flatten()))
        t_mean = summed_norm_cnts.mean()
        t_sigma = loess_gene.findy(t_mean)
        gene = test_genes[ig] 
        #test every testable SNVs and write output 
        #extract SNVs
        snv_coords = data_testable["keys"][t_index]
        snv_ratios = ars #data_testable['allelic_ratio'][t_index]
        REF_snv_ratios = REF_ars
        snv_regions= arrayed_test_matrix[t_index][:,ig] 
        snv_sigma  = data_testable['expected_sigma'][t_index]
        snv_major  = data_testable['XM'][t_index]
        snv_minor  = data_testable['Xm'][t_index]
        snv_padjs  = data_testable['adjusted_P'][t_index]
        snv_porgs  = data_testable['original_P'][t_index]
        assert len(snv_coords) == len(snv_ratios)
        assert len(snv_coords) == len(snv_regions)
        assert len(snv_coords) == len(snv_sigma)
        assert len(snv_coords) == len(snv_major)
        assert len(snv_coords) == len(snv_minor)
        assert len(snv_coords) == len(snv_padjs)
        assert len(snv_coords) == len(snv_porgs)
        assert t_norm_cnts.shape == t_raw_cnts.shape
        assert t_norm_cnts.shape[0] == len(snv_coords)
        #ouput for gene part:
        # Gene, #testable SNV,  #significant SNV, #SNV<Nsigma,  mean coverage, expected sigma, allelic_Ratio (ASE), allelic_ratio (testable), information of sites
        npow = len(t_index)
        nsig_old = len(psig_index)
        nsig_new = len(np.where(np.abs(delta/t_sigma) < opts.nsigma)[0])
        cov = t_mean
        sigma_exp = t_sigma
        #output for SNV part:
        #Site, coverage (R1+R2+A1+A2/4), ASE (pvalue or adjusted pvalue), expected sigma, delta_ar,[coverage(R1/A1, R2/A2), normed coverage(R1/A1, R2/A2)]; 
        snv_vals =[]
        for i in range(len(snv_coords)):
            crd = snv_coords[i]
            ar  = snv_ratios[i]
            REF_ar = REF_snv_ratios[i]
            dar = ar-avg_ar
            ars = "{:.2f}".format(ar)
            REF_ars = "{:.2f}".format(REF_ar)
            dars = "{:.2f}".format(dar)
            reg = PRIORITY_CODES[snv_regions[i]]
            sigma = "{:.2e}".format(snv_sigma[i])
            major = "{:.2f}".format(snv_major[i])
            minor = "{:.2f}".format(snv_minor[i])
            porg= "{:.2e}".format(snv_porgs[i])
            padj= "{:.2e}".format(snv_padjs[i])
            nraw = "/".join([str(i) for i in t_raw_cnts[i].flatten()])
            nnorm= "/".join(["{:.2f}".format(i) for i in t_norm_cnts[i].flatten()])
            val=str(crd)+","+str(reg)+","+str(ars)+"," + str(REF_ars) + ","  +str(dars)+","+str(sigma)+","+str(major)+","+str(minor)+","+str(porg)+","+str(padj)+","+nraw+","+nnorm
            snv_vals.append(val)
        assert len(snv_vals) == npow
        if np.all(np.abs(delta/t_sigma) < opts.nsigma):
            foutput.write("{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.3f}\t{:.3f}\t{:.3f}\t".format (gene, "ASE",    npow, nsig_old, nsig_new, cov, sigma_exp, avg_ar_sig, avg_ar))
            for snv_val in snv_vals:
                foutput.write("{};".format(snv_val))
            foutput.write("\n")
            N_ASE += 1
        else:
            foutput.write("{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.3f}\t{:.3f}\t{:.3f}\t".format (gene, "nonASE", npow, nsig_old, nsig_new, cov, sigma_exp, avg_ar_sig, avg_ar))
            for snv_val in snv_vals:
                foutput.write("{};".format(snv_val))
            foutput.write("\n")
        bar.next() 
logger.critical("Number of ASE genes is {}".format(N_ASE))
#print("Number of ASE genes is {}".format(N_ASE))

arrayed_powerful_matrix=powerful_matrix.toarray()
logger.info("I am working on Loess with powerful genes, N={}".format(arrayed_powerful_matrix.shape[1]))
for ig in range(arrayed_powerful_matrix.shape[1]):
    t_index = np.where(arrayed_powerful_matrix[:,ig] > 0 )[0]
    assert len(t_index) >= 2
    t_norm_cnts = data_testable['normalized_counts'][t_index]
    t_raw_cnts = data_testable['counts'][t_index]
    ##way1 for average allelic ratio of the gene: using only the ASE SNVs
    ##extract adjusted P values of each SNV in the gene and then extract the sites with significant adjusted P values
    psig_index = np.where(data_testable["adjusted_P"][t_index]<opts.p)[0]
    assert float(len(psig_index))<=0.5*float(len(t_index))
    sig_norm_cnts =t_norm_cnts[psig_index]
    if len(psig_index)>0:
        summed_sig_norm_cnts = sig_norm_cnts.sum(axis=1)
        ars_sig = summed_sig_norm_cnts.max(axis=1) / summed_sig_norm_cnts.sum(axis=1)
        avg_ar_sig = np.mean(ars_sig)
    else:
        avg_ar_sig = "NA"
    #way2 for average allelic ratio of the gene: using only all the testable SNVs
    #calculate delta allelic ratio of all testable SNVs in the gene
    summed_norm_cnts = t_norm_cnts.sum(axis=1) #shape nsites (in the gene) X 2 , sum the replicates
    ars = summed_norm_cnts.max(axis=1) / summed_norm_cnts.sum(axis=1)
    REF_ars = summed_norm_cnts[:,0] / summed_norm_cnts.sum(axis=1)
    avg_ar = np.mean(ars)
    delta = ars - avg_ar
    t_mean = summed_norm_cnts.mean()
    gene = powerful_genes[ig] 
    #write every powerful genes to output 
    #extract SNVs
    snv_coords = data_testable["keys"][t_index]
    snv_ratios = ars 
    REF_snv_ratios = REF_ars
    snv_deltas = delta
    snv_regions= arrayed_powerful_matrix[t_index][:,ig] 
    snv_sigma  = data_testable['expected_sigma'][t_index]
    snv_major  = data_testable['XM'][t_index]
    snv_minor  = data_testable['Xm'][t_index]
    snv_padjs  = data_testable['adjusted_P'][t_index]
    snv_porgs  = data_testable['original_P'][t_index]
    assert len(snv_coords) == len(snv_ratios)
    assert len(snv_coords) == len(snv_regions)
    assert len(snv_coords) == len(snv_sigma)
    assert len(snv_coords) == len(snv_major)
    assert len(snv_coords) == len(snv_minor)
    assert len(snv_coords) == len(snv_padjs)
    assert len(snv_coords) == len(snv_porgs)
    assert t_norm_cnts.shape == t_raw_cnts.shape
    assert t_norm_cnts.shape[0] == len(snv_coords)

    #ouput for gene part:
    # Gene, #testable SNV,  #significant SNV, #SNV<Nsigma,  mean coverage, expected sigma, allelic_Ratio (ASE), allelic_ratio (testable), information of sites
    npow = len(t_index)
    nsig_old = len(psig_index)
    nsig_new = "NA"
    cov = t_mean
    sigma_exp = "NA" #t_sigma
    #output for SNV part:
    snv_vals =[]
    for i in range(len(snv_coords)):
        crd = snv_coords[i]
        ar  = snv_ratios[i]
        dar = snv_deltas[i]
        REF_ar = REF_snv_ratios[i]
        ars = "{:.2f}".format(ar)
        REF_ars = "{:.2f}".format(REF_ar)
        dars = "{:.2f}".format(dar)
        reg = PRIORITY_CODES[snv_regions[i]]
        sigma = "{:.2e}".format(snv_sigma[i])
        major = "{:.2f}".format(snv_major[i])
        minor = "{:.2f}".format(snv_minor[i])
        porg= "{:.2e}".format(snv_porgs[i])
        padj= "{:.2e}".format(snv_padjs[i])
        nraw = "/".join([str(i) for i in t_raw_cnts[i].flatten()])
        nnorm= "/".join(["{:.2f}".format(i) for i in t_norm_cnts[i].flatten()])
        val=str(crd)+","+str(reg)+","+str(ars)+","+str(REF_ars)+","+str(dars)+","+str(sigma)+","+str(major)+","+str(minor)+","+str(porg)+","+str(padj)+","+nraw+","+nnorm
        snv_vals.append(val)
    assert len(snv_vals) == npow

    #Site, coverage (R1+R2+A1+A2/4), ASE (pvalue or adjusted pvalue), expected sigma, delta_ar,[cov(R1/A1, R2/A2), normed cov(R1/A1, R2/A2)]; 
    foutput.write("{}\t{}\t{}\t{}\t{}\t{:.1f}\t{}\t{}\t{:.3f}\t".format (gene, "nonASE", npow, nsig_old, nsig_new, cov, sigma_exp, avg_ar_sig, avg_ar))
    for snv_val in snv_vals:
        foutput.write("{};".format(snv_val))
    foutput.write("\n")
foutput.close()

if opts.plotsigma:
    #check the distribution of sigma (of delta allelic ratio) in control and testable genes
    logger.critical("Number of control genes: {}".format(len(control_genes)))
    logger.critical("Number of testable genes: {}".format(len(test_genes)))
    logger.critical('number of ASE genes {}'.format(len(ase_genes)))
    logger.critical('number of ASE-part genes {}'.format(len(partase_genes)))
    
    #calculate the mean and sigma of delta allelic ratios of (previouse defined) ASE genes and part-ASE genes
    #mean allelic raito of a gene is by all testable SNVs on the gene (not just ASE sites on the gene)
    arrayed_ase_matrix=ase_matrix.toarray()
    for ig in range(arrayed_ase_matrix.shape[1]):
        t_index = np.where(arrayed_ase_matrix[:,ig] > 0 )[0]
        t_norm_cnts = data_testable['normalized_counts'][t_index]
        #calculate delta allelic ratio of all testable SNVs in the gene (here all testable SNVs are ASE SNVs)
        summed_norm_cnts = t_norm_cnts.sum(axis=1) #shape nsites (in the gene) X 2 , sum the replicates
        ars = summed_norm_cnts.max(axis=1) / summed_norm_cnts.sum(axis=1)
        avg_ar = np.mean(ars)
        delta = ars - avg_ar
        plotting_data['ase_set_mean'].append(delta.mean())
        t0 = np.asarray([delta, -1.0 * delta])
        t = t0.flatten()
        sigma = np.std(t)
        plotting_data['ase_set_sigma'].append(sigma)

    arrayed_partase_matrix=partase_matrix.toarray()
    for ig in range(arrayed_partase_matrix.shape[1]):
        t_index = np.where(arrayed_partase_matrix[:,ig] > 0 )[0]
        t_norm_cnts = data_testable['normalized_counts'][t_index]

        #calculate the mean allelic ratio of the gene using all testable SNVs (not just ASE SNVs)
        summed_norm_cnts = t_norm_cnts.sum(axis=1) #shape nsites (in the gene) X 2 , sum the replicates
        ars = summed_norm_cnts.max(axis=1) / summed_norm_cnts.sum(axis=1)
        avg_ar = np.mean(ars)
        delta = ars - avg_ar
        plotting_data['partase_set_mean'].append(delta.mean())
        t0 = np.asarray([delta, -1.0 * delta])
        t = t0.flatten()
        sigma = np.std(t)
        plotting_data['partase_set_sigma'].append(sigma)

    import matplotlib.pyplot as plt
    plt.clf()
    for k,v in plotting_data.items():
        if "sigma" in k:
            bins = np.linspace(0.0,0.5,60)
            cnts = np.histogram(v,bins=bins, density=True)[0]
            plt.plot(bins[:-1],cnts,label=k)
            plt.legend()
            plt.savefig("{}_sigma_distribution.png".format(outall_prefix))
    plt.clf()
 