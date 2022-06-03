
import sys
import matplotlib.pyplot as plt
from statistics import mean,median,stdev

from functions import * # imports from functions.py (must be in the same directory as current script)

### Inputs
vcf_file = '/Disk/user/mutect2.vcf'
GTF_file = '/Disk/user/reference/dmel-all-chromosome-r6.12.gtf'
masking_file = '/Disk/user/reference/dmel-all-chromosome-r6.12.fasta.out'
###/Inputs

### Outputs
output_mainDir = '/Disk/user/mutect2_analysis_output'
###/Outputs

### GLOBALS
global_vars_modules.global_idx_iterator = 0 # keeps track of added index via assignIdx function
chroms = ['2L','2R','3L','3R','X','4']
autosomes = ['2L','2R','3L','3R']
###/GLOBALS


def plot_keyVal_arr(inp_arr,title='',xticksrotate_deg=0,ylim=[],yscale=False,savefig='',selectKey=None,separateTimelines=False):
    if type(inp_arr) == dict:       inp_arr = inp_arr.items()
    
    yvals = []
    xvals = []
    for key,val in inp_arr:
        if selectKey!=None:     val = val[selectKey]
        yvals.append(val)
        xvals.append(key)
    
    fig = plt.figure(figsize=(11,11))
    ax = plt.gca()
    bars = ax.bar(xvals,yvals)
    if title:       plt.title(title)
    if xticksrotate_deg:        ax.set_xticklabels(xvals,rotation=xticksrotate_deg)
    plt.tight_layout()
    if ylim:        plt.ylim(ylim)
    if yscale:      plt.yscale(yscale)
    
    if separateTimelines:
        trans = ax.get_xaxis_transform()
        ax2 = plt.twinx()
        ax2.set_xlim(ax.get_xlim())
        ax2.get_yaxis().set_visible(False)
        for enum,bar in enumerate(bars):
            if enum == 0: continue #backwards checking
            cur_annot = xvals[enum]
            prev_annot = xvals[enum-1]
            # check if timeline ID chaged, then paint
            if '-'.join(cur_annot.split('-')[0:2]) != '-'.join(prev_annot.split('-')[0:2]):
                barX = bar.get_x()
                width = bar.get_width()
                ax2.plot([barX-width*0.12,barX-width*0.12],[-0.1,1], linewidth=1, color='red',transform=trans,clip_on=False)
            #/
        
    if savefig:     plt.savefig(savefig)
    
######## SCRIPT START

SV_idx_map = {}
SV_len_thresh = 0
print('Import VCF file...')
vcf_fo = pyvcf.Reader(filename=vcf_file)
for entry in vcf_fo:
    # Format of VCF:
    # chrom1,pos1 -> [chrom2+pos2,...,chromN+posN]
    
    #if entry.INFO['SVTYPE'] == 'DUP': sys.exit()
    chrom1 = entry.CHROM
    pos1 = entry.POS
    pos2 = entry.end
    if pos2 == pos1: pos2 += 1
    
    ## Parse genotypes, by each ALT and REF
    gts = entry.aaf
    DP = entry.INFO['DP']
    
    # Skip if this position has support in embryo
    if entry.samples[-1].sample != 'emb': sys.exit('no emb data!!')
    if sum(entry.samples[-1].data.AD[1:]) > 0: continue
    #/
    
    # Parse supp per ALT, per sample
    sample_supps = {}
    alt_seqs = {}
    for enum,ALT in enumerate(entry.ALT):
        alt_seqs[enum] = ALT.sequence
        for chunk in entry.samples:
            sample = chunk.sample
            data = chunk.data
            
            init(sample_supps,sample,{'AOs':[],'AFs':[]})
            sample_supps[sample]['AOs'].append(data.AD[enum+1]) # 1 is "ref supp"
            if type(data.AF) == list:
                sample_supps[sample]['AFs'].append(data.AF[enum])
            else:
                sample_supps[sample]['AFs'].append(data.AF)
            sample_supps[sample]['RO'] = data.AD[0]
            
            try:
                sample_supps[sample]['DP'] = data.DP
            except:
                print('Failed to grab DP! Using DP as entry.INFO[DP]')
                sample_supps[sample]['DP'] = DP
    #/
    
    tmp_SV = {'rname':chrom1,'rcoords':[pos1,pos2],'ref':entry.REF,'alts':alt_seqs,'supps':sample_supps}
    assignIdx(tmp_SV)
    SV_idx_map[tmp_SV['idx']] = tmp_SV
###/ Import VCFs

### Parse basic information
## Get all sample names
SAMPLES = {}
for SV_idx,SV in SV_idx_map.items():
    for TLid in SV['supps']:
        init(SAMPLES,TLid,0)
        SAMPLES[TLid] += 1
SAMPLES_noEmb = set(SAMPLES).difference({'emb'})
##/
## Parse timeline samples
TIMELINES_samples = {}
SAMPLES_timeline = {}
for sample in SAMPLES:
    timeline = '-'.join(sample.split('-')[:2])
    init(TIMELINES_samples,timeline,set())
    TIMELINES_samples[timeline].add(sample)
    SAMPLES_timeline[sample] = timeline
##/
## Parse timeline samples by age
SAMPLES_age = {}
AGE_samples = {}
for sample in SAMPLES:
    if sample == 'emb': continue
    time = int(sample.split('-')[-1])
    age = None
    if getRangeOvlp([time]*2,[0,150]) >= 0:      age = '1'
    elif getRangeOvlp([time]*2,[150,300]) >= 0:  age = '2'
    elif time > 300:                             age = '3'
    if not age: sys.exit('unable to assign age. check me!')
    SAMPLES_age[sample] = age
    
    init(AGE_samples,age,set())
    AGE_samples[age].add(sample)
##/
    
## Roughly parse depth per sample
SAMPLES_rname_depth_pre = {}
for SV_idx,SV in SV_idx_map.items():
    rname = SV['rname']
    if not rname in chroms: continue
    for sample,data in SV['supps'].items():
        DP = data['DP']
        if getRangeOvlp([DP]*2,[10,100]) >= 0:
            init(SAMPLES_rname_depth_pre,sample,{})
            init(SAMPLES_rname_depth_pre[sample],rname,[])
            SAMPLES_rname_depth_pre[sample][rname].append(DP)

SAMPLES_rname_depth = {}
for sample,rname_depths in SAMPLES_rname_depth_pre.items():
    for rname,depths in rname_depths.items():
        init(SAMPLES_rname_depth,sample,{})
        SAMPLES_rname_depth[sample][rname] = median(depths)
##/
## Depth and X to A ratio
    SAMPLES_X_TO_A_ratios = {}
    SAMPLES_AVGCOV = {}
    for sample,rname_depths in SAMPLES_rname_depth_pre.items():
        A = mean(rname_depths['2L'] + rname_depths['2R'] + rname_depths['3L'] + rname_depths['3R'])
        X = mean(rname_depths['X'])
        X_to_A = X/A
        
        SAMPLES_X_TO_A_ratios[sample] = X_to_A
        SAMPLES_AVGCOV[sample] = A
##/
###/

### Add masking + gtf overlaps to SVs
# Import masking
masking = importRepeatMasker(masking_file,rSorted=True)
masking_idx_map = {}
for rname in masking:
    for entry in masking[rname]:
        assignIdx(entry)
        masking_idx_map[entry['idx']] = entry
masking_idx_map_ROMD = rangeOverlaps_makeMappingDict(masking_idx_map,100,coordsKey='rcoords',sortByKey='rname')

# Import GTFs
GTF_data = importGTF(GTF_file,rSorted=True)
GTF_types = set()
GTF_idx_map = {}
gene_GTF_idx_map = {}
for rname in GTF_data:
    for entry in GTF_data[rname]:
        assignIdx(entry)
        GTF_idx_map[entry['idx']] = entry
        GTF_types.add(entry['type'])
        
        if entry['type'] == 'gene':
            gene_GTF_idx_map[entry['gene_id']] = entry['idx']
GTF_idx_map_ROMD = rangeOverlaps_makeMappingDict(GTF_idx_map,1000,coordsKey='rcoords',sortByKey='rname')

# Sort GTFs by gene, sort GTFs by type
# gene -> GTFs, GTF -> gene
gene_GTFs = {}
GTFs_gene = {}
for gtf_idx,gtf_entry in GTF_idx_map.items():
    gene = gtf_entry['gene_id']
    
    # Gene -> GTF idxs
    init(gene_GTFs,gene,[])
    gene_GTFs[gene].append(gtf_idx)
    
    # GTFs -> Gene
    for entry in gtf_entry['data'].split(';'):
        if entry.find('FB') != -1:
            entry2 = entry.split()
            entry2[1] = entry2[1].replace('"','')
            init(GTFs_gene,entry2[1],gene)
#/

## Add masking + gtf to SVs
for SV_idx,SV in SV_idx_map.items():
    rname = SV['rname']
    rcoords = SV['rcoords']
    
    # Check mask overlaps
    maskClasses = set()
    ovlps = rangeOverlaps_lookup(rcoords,masking_idx_map_ROMD,100,map_lookup_key=rname)
    maskRanges = []
    for oS,mask_idx,oE in ovlps:
        maskE = masking_idx_map[mask_idx]
        typee = maskE['class']
        maskClasses.add(typee)
        if oS != oE:
            maskRanges.append([oS,mask_idx,oE])
        else:
            maskRanges.append([oS,mask_idx,oS+1])
    ovlpRanges = computeOvlpRanges_wIdxs(maskRanges)
    mask_numBP = 0
    for oS,idxs,oE in ovlpRanges:
        mask_numBP += oE-oS
    mask_covFrac = mask_numBP / (rcoords[-1]-rcoords[0])
    #/
    
    # Check GTF overlaps
    gtf_ovlps = {}
    gtfClasses = set()
    ovlps = rangeOverlaps_lookup(rcoords,GTF_idx_map_ROMD,1000,map_lookup_key=rname)
    for oS,gtf_idx,oE in ovlps:
        GTF = GTF_idx_map[gtf_idx]
        typee = GTF['type']
        gene_id = GTF['gene_id']
        
        init(gtf_ovlps,gene_id,{})
        init(gtf_ovlps[gene_id],typee,set())
        gtf_ovlps[gene_id][typee].add(gtf_idx)
        
        gtfClasses.add(typee)
    #/
    
    # Save
    SV['masked_covFrac'] = mask_covFrac
    SV['maskClasses'] = maskClasses
    SV['gtf_ovlps'] = gtf_ovlps
    SV['gtfClasses'] = gtfClasses
    #/
##/
###/
    
### Summarize sample datasets: Fraction of dataset with X% AF
# find fraction of each dataset that have SNPs in e.g. AF 0.6-0.7
samples_AFs_fractions = {}
AF_step = 0.1

# Init deviate bins
for any_SV in SV_idx_map.values():
    for sample in any_SV['supps']:
        init(samples_AFs_fractions,sample,{})
        val = 0
        while val <= 1:
            val = round(val,2)
            init(samples_AFs_fractions[sample],val,0) #round to remove bullshit 1.1999999999999 instead of 1.2
            val += AF_step
    break
#/

# Clear AF bin marker key
SK = 'sample_AF_bins'
for SV in SV_idx_map.values():
    if SK in SV:        del SV[SK]
#/

for SV in SV_idx_map.values():
    # Run filters
    if 1 and 'skip SV if (1) not single-nucleotide':
        SVlens = []
        for altenum,alt in SV['alts'].items():
            SVlen = abs(len(alt)-len(SV['ref']))
            SVlens.append(SVlen)
        if max(SVlens) >= 2: continue #skip if not a SNP
        
    if 1 and 'skip if masked':
        if SV['masked_covFrac'] > 0: continue
    
    if 1 and 'skip if below X DP':
        samples_below_DP = 0
        for sample,supp in SV['supps'].items():
            if supp['DP'] < 10:
                samples_below_DP += 1
        
        if samples_below_DP > 2: continue # skip if more than 2 samples had less than 10 DP here...
    #/
    # Add data
    for sample,supp in SV['supps'].items():
        AO_sum = sum(supp['AOs'])
        if AO_sum == 0: continue #skip if no read support for SNP
        
        AF_sum = sum(supp['AFs'])
        
        AF_bin = round(round(AF_sum/AF_step)*AF_step,2) # outer round is to remove bullshit 1.199999999 instead of 1.2
        
        if AF_bin == 0: continue #skip 0-bin
        
        # Save marker at SV
        SK = 'sample_AF_bins'
        init(SV,SK,{})
        SV[SK][sample] = AF_bin
        #/
        
        # save
        samples_AFs_fractions[sample][AF_bin] += 1
        #/
    #/

# Transform absolute numbers to density
samples_AFs_fractions_proportion = {}
for sample in samples_AFs_fractions:
    for deviate_bin in samples_AFs_fractions[sample]:
        init(samples_AFs_fractions_proportion,sample,{})
        samples_AFs_fractions_proportion[sample][deviate_bin] = round(samples_AFs_fractions[sample][deviate_bin] / sum(list(samples_AFs_fractions[sample].values())),6)
#/

# Dump?
if 0 and 'output for excel':
    with open(output_mainDir+'/'+'SNPs_density.tsv','w') as nf:
        # Write header
        writeArr = ['sample']
        for any_sample in samples_AFs_fractions_proportion:
            for deviate_bin in sorted(samples_AFs_fractions_proportion[any_sample]):
                writeArr.append(deviate_bin)
            break # only take one sample
        writeArr.append('num') # total observations per sample
        writeArr.append('depth') # estimate of read depth per sample
        nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
        # Write sample data
        for sample in sorted(samples_AFs_fractions_proportion):
            writeArr = [sample]
            for deviate_bin in sorted(samples_AFs_fractions_proportion[sample]):
                writeArr.append(samples_AFs_fractions_proportion[sample][deviate_bin])
            writeArr.append(sum(list(samples_AFs_fractions[sample].values()))) # write total observations
            writeArr.append(SAMPLES_AVGCOV[sample]) # write depth
            nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
#/

## Compile genes per sample based on AF bin (using marked SVs from above)
# Compile genes
AF_bin_out_thresh = 0.3
samples_genes = {}
for SV in SV_idx_map.values():
    key_to_check = 'sample_AF_bins'
    if not key_to_check in SV: continue # skip if marker not present at SV
    
    for sample,AF_bin in SV[key_to_check].items():
        if AF_bin >= AF_bin_out_thresh:
            for gene,data in SV['gtf_ovlps'].items():
                
                init(samples_genes,sample,set())
                samples_genes[sample].add(gene)
        
#/
# Dump?
if 0 and 'output found genes for SVs':
    with open(output_mainDir+'/'+'SNPs_AFbin_genes.tsv','w') as nf:
        # write header
        writeArr = ['sample','genes: transpose in excel']
        nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
        # Write sample and genes
        for sample,genes in sorted(samples_genes.items(),key=lambda x: x[0]):
            nf.write(sample+'\t'+'\t'.join(sorted(genes))+'\n')
        #/
#/
##/

## Compile timeline-AF-trends between samples based on AF bin (using marked SVs for above)
# Compile SVs per timeline where any sample has >0.3 AF for the SNP
AF_bin_thresh = 0.3
SNPs_per_timeline = {}
for SV_idx,SV in SV_idx_map.items():
    key_to_check = 'sample_AF_bins'
    if not key_to_check in SV: continue # skip if marker not present at SV
    
    for sample,AF_bin in SV[key_to_check].items():
        if AF_bin >= AF_bin_thresh:
            timeline = SAMPLES_timeline[sample]
            
            # Skip timeline if it does not have at least 3 samples
            if not len(TIMELINES_samples[timeline]) >= 3: continue
            #/
            
            # Check so all timelines actually have read support for SNP (Mutect2 can yield AF<>0 with AO=0)
            samples_AOs = []
            for sample2 in TIMELINES_samples[timeline]:
                samples_AOs.append(sum(SV['supps'][sample2]['AOs']))
            if min(samples_AOs) == 0: continue #skip if not all samples has read support at SNP
            #/
            
            init(SNPs_per_timeline,timeline,{})
            init(SNPs_per_timeline[timeline],SV_idx,{})
            
            for sample2 in TIMELINES_samples[timeline]:
                SNPs_per_timeline[timeline][SV_idx][sample2] = sum(SV['supps'][sample2]['AFs'])
# plot?
if 0 and 'plot?':
    genes_updown_classi = {'up':set(),'down':set()}
    for timeline,SVs_data in SNPs_per_timeline.items():
        # construct X-axis
        xaxis = []
        for sample in sorted(TIMELINES_samples[timeline], key=lambda x: int(x.split('-')[-1])):
            xaxis.append(sample)
        #/
        # Compile all lines to plot for timeline
        plotArrs = []
        SV_idx_trace = []
        for SV_idx,data in SVs_data.items():
            tmp_arr = []
            for sample in xaxis:
                tmp_arr.append(data[sample])
            plotArrs.append(tmp_arr)
            SV_idx_trace.append(SV_idx)
        #/
        # Plot timeline
        fig = plt.figure(figsize=(11,11))
        ax = plt.gca()
        num_increase = 0
        num_decrease = 0
        for enum,plotArr in enumerate(plotArrs):
            # for genes_updown_classi: get SV gene ovlp
            SV_genes = set()
            for gene,ovlptypes in SV_idx_map[SV_idx_trace[enum]]['gtf_ovlps'].items():
                if set(ovlptypes).intersection(set(['gene','mRNA'])) == set(ovlptypes): continue #skip gene if it ovlp in intron
                SV_genes.add(gene)
            #/
            
            # set color based on increase over time or decrease over time
            color = 'black'
            if plotArr[0] < plotArr[-1]:
                color = 'green'
                num_increase += 1
                genes_updown_classi['up'].update(SV_genes)
            elif plotArr[0] >= plotArr[-1]:
                color = 'red'
                num_decrease += 1
                genes_updown_classi['down'].update(SV_genes)
            ax.plot(plotArr,color=color)
            
        ax.set_xticklabels(xaxis)
        ax.set_xticks(list(range(len(xaxis))))
        plt.title(timeline + '||' + '+'+str(num_increase) + ' -'+str(num_decrease))
        plt.savefig(output_mainDir+'/'+'SNPs_overtimelines.'+timeline+'.pdf')
        #/
        
    if 0 and 'dump up/down gene lists?':
        with open(output_mainDir+'/'+'genes_up_overtimelines.tsv','w') as nf:
            nf.write('\n'.join(sorted(genes_updown_classi['up'])))
        with open(output_mainDir+'/'+'genes_down_overtimelines.tsv','w') as nf:
            nf.write('\n'.join(sorted(genes_updown_classi['down'])))
#/
##/
###/