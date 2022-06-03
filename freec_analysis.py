
import sys
import os
import matplotlib.pyplot as plt
from statistics import mean,median,stdev

from functions import * # imports from functions.py (must be in the same directory as current script)

### Inputs
masking_file = '/Disk/user/reference/dmel-all-chromosome-r6.12.fasta.out'
GTF_file = '/Disk/user/reference/dmel-all-chromosome-r6.12.gtf'

samples_dir = '/Disk/user/FREEC' # directory to scan for folders of sample FREEC outputs
freec_basename = '<sample>/freec_default/alignments.sorted.bam_ratio.BedGraph'  # will replace <sample> with folder names (of samples) in "samples_dir"
###/

### OUTPUT
output_dir = '/Disk/user/FREEC_analysis_output'
###/

### Define globals
chroms = ['2L','2R','3L','3R','X']

global_vars_modules.global_idx_iterator = 0 # keeps track of added index via assignIdx function

############ SCRIPT START

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
        
### Import Ref stuff
## Masking
masking = importRepeatMasker(masking_file,rSorted=True)
masking_idx_map = {}
for rname in masking:
    for entry in masking[rname]:
        assignIdx(entry)
        masking_idx_map[entry['idx']] = entry
masking_idx_map_ROMD = rangeOverlaps_makeMappingDict(masking_idx_map,100,coordsKey='rcoords',sortByKey='rname')

def calc_ovlp_numBP_frac(rname,rcoords,ROMD,ROMD_step):
    ovlps = rangeOverlaps_lookup(rcoords,ROMD,ROMD_step,map_lookup_key=rname)
    maskRanges = []
    for oS,mask_idx,oE in ovlps:
        if oS != oE:
            maskRanges.append([oS,mask_idx,oE])
        else:
            maskRanges.append([oS,mask_idx,oS+1])
    ovlpRanges = computeOvlpRanges_wIdxs(maskRanges)
    mask_numBP = 0
    for oS,idxs,oE in ovlpRanges:
        mask_numBP += oE-oS
    mask_covFrac = mask_numBP / (rcoords[-1]-rcoords[0])
    
    return mask_numBP,mask_covFrac
##/
## Import GTFs
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
##/
###/

### Import FREEC & assign masking
rname_freecs = {}
for sample in os.listdir(samples_dir):
    freec_file = samples_dir+'/'+freec_basename.replace('<sample>',sample)
    if not os.path.exists(freec_file): continue

    bin_size = 1000
    with open(freec_file,'r') as f:
        for ln,line in enumerate(f): #Use ln to assign IDs
            line = line.strip('\n')
            line = line.split()
            if line[0] == 'track': continue
            #rname,start,ratio,ratio_median,CN = line
            rname,start,end,ratio = line
            rname = rname.replace('chr','') #remove added chromosome identifier
            if not rname in ('2L','2R','3L','3R','X','4','Y','mitochondrion_genome'): continue
            rstart=int(start)
            rend=int(end)-1
            ratio=float(ratio)
            
            # Init freec if not existant (freec-entry will hold every sample value)
            init(rname_freecs,rname,{})
            if not rstart in rname_freecs[rname]:
                rcoords = [rstart,rend]
                maskNumBP,maskCovFrac = calc_ovlp_numBP_frac(rname,rcoords,masking_idx_map_ROMD,100)
                
                tmp_freecE = {'rname':rname,'rcoords':rcoords,'samples':{},'masked_covFrac':maskCovFrac}
                rname_freecs[rname][rstart] = tmp_freecE
            #/
            
            # Update existing freec with sample data
            rname_freecs[rname][rstart]['samples'][sample] = ratio
            #/
###/

### Parse samples & timelines & age
## Parse samples
SAMPLES = {}
for rname,freecs in rname_freecs.items():
    rname_sample_data = {}
    for _,freec in freecs.items():
        for sample in freec['samples']:
            init(SAMPLES,sample,0)
            SAMPLES[sample] += 1
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
###/
            
### INFO: Plot chromosome means/medians
if 0:
    rname_sample_stats = {}
    for rname,freecs in rname_freecs.items():
        rname_sample_data = {}
        for _,freec in freecs.items():
            if freec['masked_covFrac'] >= 0.3: continue
            for sample,ratio in freec['samples'].items():
                init(rname_sample_data,sample,[])
                rname_sample_data[sample].append(ratio)
        
        init(rname_sample_stats,rname,{})
        for sample,ratios in rname_sample_data.items():
            rname_sample_stats[rname][sample] = {'mean':mean(ratios),'median':median(ratios)}
    
    fig = plt.figure(figsize=(11,11))
    xlabels = []
    yvals = []
    emb_vals = []
    for rname,sample_data in sorted(rname_sample_stats.items(),key=lambda x: x[0]):
        if rname == 'mitochondrion_genome': continue
        rname_data = []
        for sample,data in sample_data.items():
            if not sample == 'emb':
                rname_data.append(data['mean'])
            else:
                emb_vals.append(data['mean'])
        xlabels.append(rname)
        yvals.append(rname_data)
    
    """
    for enum,_ in enumerate(yvals):
        plt.scatter([enum]*len(yvals[enum]),yvals[enum])
    """
    
    plt.boxplot(yvals)
    plt.scatter([1,2,3,4,5,6,7],emb_vals,color='red')
    ax = plt.gca()
    #ax.set_xticks(list(range(len(xlabels))))
    ax.set_xticklabels(xlabels)
    plt.show()
###/
        
### For each sample, make chromosomal arrays of freec assignment
repma_cutoff = 0.1 #skip freec bins with more than this fraction overlapping repetitive sequence
sample_rname_blocks = {}
for rname in rname_freecs:
    if rname == 'mitochondrion_genome': continue
    for rstart,freec in rname_freecs[rname].items():
        for sample in freec['samples']:
            # Skip if repeatmasked
            if freec['masked_covFrac'] >= repma_cutoff: continue
            #/
            init(sample_rname_blocks,sample,{})
            init(sample_rname_blocks[sample],rname,[])
            sample_rname_blocks[sample][rname].append({'rcoords':freec['rcoords'],'rname':rname,'sample':sample,
                                                       'ratio':freec['samples'][sample],'masked_covFrac':freec['masked_covFrac']})

# Sort blocks by coord
for sample,rname_blocks in sample_rname_blocks.items():
    for rname,blocks in rname_blocks.items():
        blocks = sorted(blocks,key=lambda x: x['rcoords'][0])
#/

# Find sample chrom baseline stats (e.g. mean, median, stdev)
sample_rname_stats = {}
for sample in sample_rname_blocks:
    for rname in sample_rname_blocks[sample]:
        block_ratios = []
        for block in sample_rname_blocks[sample][rname]:
            if block['ratio'] > 10: continue #skip if above 10 in copy number
            block_ratios.append(block['ratio'])
            
        if len(block_ratios) < 2: continue #skip if not at least two blocks persisted
        
        blocks_mean = mean(block_ratios)
        blocks_med = median(block_ratios)
        blocks_stdev = stdev(block_ratios)
        
        init(sample_rname_stats,sample,{})
        init(sample_rname_stats[sample],rname,{})
        sample_rname_stats[sample][rname]['mean'] = blocks_mean
        sample_rname_stats[sample][rname]['med'] = blocks_med
        sample_rname_stats[sample][rname]['stdev'] = blocks_stdev
#/
###/
        
### Task1: Summarize sample datasets: Fraction of dataset deviating X% from baseline
samples_CNratio_fractions = {}
deviate_max = 2
deviate_step = 0.2

# Init deviate bins
for sample,rname_blocks in sample_rname_blocks.items():
    init(samples_CNratio_fractions,sample,{})
    val = -deviate_max
    while val <= deviate_max:
        val = round(val,2)
        init(samples_CNratio_fractions[sample],val,0) #round to remove bullshit 1.1999999999999 instead of 1.2
        val += deviate_step
#/

# Fill deviate bins
for sample,rname_blocks in sample_rname_blocks.items():
    for rname,blocks in rname_blocks.items():
        emb_baseline = sample_rname_stats['emb'][rname]['mean']
        for block in blocks:
            block_deviate = block['ratio'] - emb_baseline
            
            if abs(block_deviate) >= deviate_max: continue #skip if very high ratio
        
            # Calc deviate-bin
            deviate_bin = round(round(block_deviate/deviate_step)*deviate_step,2) # outer round is to remove bullshit 1.199999999 instead of 1.2
            if deviate_bin > deviate_max:          deviate_bin = deviate_max
            elif deviate_bin < -deviate_max:       deviate_bin = -deviate_max
            #/
            
            # Save count
            samples_CNratio_fractions[sample][deviate_bin] += 1
            #/
#/

# Transform absolute numbers to density
samples_CNratio_fractions_proportion = {}
for sample in samples_CNratio_fractions:
    for deviate_bin in samples_CNratio_fractions[sample]:
        init(samples_CNratio_fractions_proportion,sample,{})
        samples_CNratio_fractions_proportion[sample][deviate_bin] = round(samples_CNratio_fractions[sample][deviate_bin] / sum(list(samples_CNratio_fractions[sample].values())),6)
#/

# Dump?
if 0 and 'output for excel':
    with open(output_dir+'/'+'blocks_density.tsv','w') as nf:
        # Write header
        writeArr = ['sample']
        for any_sample in samples_CNratio_fractions_proportion:
            for deviate_bin in sorted(samples_CNratio_fractions_proportion[any_sample]):
                writeArr.append(deviate_bin)
            break # only take one sample
        nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
        # Write sample data
        for sample in sorted(samples_CNratio_fractions_proportion):
            writeArr = [sample]
            for deviate_bin in sorted(samples_CNratio_fractions_proportion[sample]):
                writeArr.append(samples_CNratio_fractions_proportion[sample][deviate_bin])
            nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
###/
    
### Task2: Find copy-number gain/losses
# make function to find GTF overlaps with blocks (can also take single blocks)
def get_blocks_gtfs(arr_of_blocks,gtfs_to_scan=set()):
    blocks_ovlps = {} # gtf_idx -> data
    
    # check if single entry input
    if type(arr_of_blocks) == dict:         arr_of_blocks = [arr_of_blocks]
    
    for block in arr_of_blocks:
        rname = block['rname']
        rcoords = block['rcoords']
        
        for oS,idx,oE in rangeOverlaps_lookup(rcoords,GTF_idx_map_ROMD,1000,map_lookup_key=rname):
            gtf = GTF_idx_map[idx]
            typee = gtf['type']
            if gtfs_to_scan and not typee in gtfs_to_scan: continue
            blocks_ovlps[idx] = {'rcoords':[oS,oE],'len':oE-oS,'type':gtf['type'],'gene_id':gtf['gene_id']}
    return blocks_ovlps

def compile_sample_rname_block_dict(input_dict,gtfs_to_scan=set()):
    genes_sampleHits = {}
    for sample in input_dict:
        for rname in input_dict[sample]:
            for window in input_dict[sample][rname]:
                gtf_ovlps = get_blocks_gtfs(window,gtfs_to_scan=gtfs_to_scan)
                for gtf_idx,data in gtf_ovlps.items():
                    gene = data['gene_id']
                    typee = data['type']
                    init(genes_sampleHits,gene,{})
                    init(genes_sampleHits[gene],sample,{})
                    init(genes_sampleHits[gene][sample],typee,0)
                    genes_sampleHits[gene][sample][typee] += 1
    return genes_sampleHits
#/

# Find sample blocks below/above baseline coverage
samples_rname_blocks_gain = {}
samples_rname_blocks_loss = {}
for sample in sample_rname_blocks:
    for rname in sample_rname_blocks[sample]:
        for block_enum,block in enumerate(sample_rname_blocks[sample][rname]):
            baseline = sample_rname_stats[sample][rname]
            
            # Get embryo block + statistic
            block_emb = sample_rname_blocks['emb'][rname][block_enum]
            baseline_emb = sample_rname_stats['emb'][rname]
            #/
            
            # check if gain (and not in emb)
            if (block['ratio'] > baseline['mean']+baseline['stdev']) and not (block_emb['ratio'] > baseline_emb['mean']+baseline_emb['stdev']):
                init(samples_rname_blocks_gain,sample,{})
                init(samples_rname_blocks_gain[sample],rname,[])
                samples_rname_blocks_gain[sample][rname].append(block)
            #/
            # check if loss (and not in emb)
            if (block['ratio'] < baseline['mean']-baseline['stdev']) and not (block_emb['ratio'] < baseline_emb['mean']-baseline_emb['stdev']):
                init(samples_rname_blocks_loss,sample,{})
                init(samples_rname_blocks_loss[sample],rname,[])
                samples_rname_blocks_loss[sample][rname].append(block)
            #/
#/
# Check gtfs under block gains/losses and compare to rand
gtfs_to_scan = {'exon', 'pre_miRNA', 'snRNA', '5UTR', 'tRNA', 'ncRNA', 'start_codon', 'miRNA', 'snoRNA', 'rRNA', 'stop_codon'}
#samples_rname_windows_gain, samples_rname_windows_gain_rand, samples_rname_windows_loss, samples_rname_windows_loss_rand
genes_sampleHits_blocks_gain = compile_sample_rname_block_dict(samples_rname_blocks_gain,gtfs_to_scan=gtfs_to_scan)
genes_sampleHits_blocks_loss = compile_sample_rname_block_dict(samples_rname_blocks_loss,gtfs_to_scan=gtfs_to_scan)
#/
# Check number of blocks found up/down for each sample
samples_numBlocks_gain = {}
samples_numBlocks_loss = {}
for sample in SAMPLES_noEmb:
    blocks_gained = 0
    for rname,blocks in samples_rname_blocks_gain[sample].items():
        blocks_gained += len(blocks)
    
    blocks_loss = 0
    for rname,blocks in samples_rname_blocks_loss[sample].items():
        blocks_loss += len(blocks)
    
    samples_numBlocks_gain[sample] = blocks_gained
    samples_numBlocks_loss[sample] = blocks_loss

if 0 and 'plot':
    plot_keyVal_arr(sorted(samples_numBlocks_gain.items(),key=lambda x: x[0]),xticksrotate_deg=90,separateTimelines=True)
    plot_keyVal_arr(sorted(samples_numBlocks_loss.items(),key=lambda x: x[0]),xticksrotate_deg=90,separateTimelines=True)
#/

### Task3: Chain blocks and find if blocks generally gain/lose coverage in larger regions than 1 block
# chain_blocks, TRY2: Run sliding window and find gained/lost blocks
def classify_window_blocks_status(blocks,sample,rname):
    baseline = sample_rname_stats[sample][rname]
    blocks_status = []
    for block in blocks:
        if (block['ratio'] > baseline['mean']+baseline['stdev']):        blocks_status.append('gain')
        elif (block['ratio'] < baseline['mean']-baseline['stdev']):      blocks_status.append('loss')
        else:                                                            blocks_status.append(None)
    
    return blocks_status

def classify_window_status(blocks_status,blocks_emb_status):    
    window_status = None
    # Check if blocks are gain
    if blocks_status.count('gain') == window_size and blocks_emb_status.count('gain') == 0:
        window_status = 'gain'
    # check if blocks are loss
    if blocks_status.count('loss') == window_size and blocks_emb_status.count('loss') == 0:
        window_status = 'loss'
        
    return window_status


print('IDE: random run in loop below takes some time, controlled via variable RUNRAND [toggle:true/false]')
RUNRAND = False
window_size = 5
samples_rname_windows_gain = {}
samples_rname_windows_gain_rand = {}
samples_rname_windows_loss = {}
samples_rname_windows_loss_rand = {}
for sample in sample_rname_blocks:
    for rname in sample_rname_blocks[sample]:
        blockStart = 0
        while (blockStart+window_size) < len(sample_rname_blocks[sample][rname]):
            ## Run real
            # get blocks for window
            blocks = []
            blocks_emb = []
            while len(blocks) < window_size:
                blocks.append(sample_rname_blocks[sample][rname][ blockStart+len(blocks) ])
                blocks_emb.append(sample_rname_blocks['emb'][rname][ blockStart+len(blocks) ])
            #/
            
            # classify gain/loss status for blocks
            blocks_status = classify_window_blocks_status(blocks,sample,rname)
            blocks_emb_status = classify_window_blocks_status(blocks_emb,'emb',rname)
            #/
            
            # Classify gain/loss status for window
            window_status = classify_window_status(blocks_status,blocks_emb_status)
            #/
            
            # Save
            if window_status == 'gain':
                init(samples_rname_windows_gain,sample,{})
                init(samples_rname_windows_gain[sample],rname,[])
                samples_rname_windows_gain[sample][rname].append(blocks)
            elif window_status == 'loss':
                init(samples_rname_windows_loss,sample,{})
                init(samples_rname_windows_loss[sample],rname,[])
                samples_rname_windows_loss[sample][rname].append(blocks)
            #/
            ##/
            ## Run rand
            if RUNRAND:
                for i in range(200): #randomize X times per real observation
                    # get blocks for window
                    blocks_rand = []
                    blocks_rand_emb = []
                    while len(blocks_rand) < window_size:
                        randBlockEnum = randint(0,len(sample_rname_blocks[sample][rname])-1)
                        blocks_rand.append(sample_rname_blocks[sample][rname][randBlockEnum])
                        blocks_rand_emb.append(sample_rname_blocks['emb'][rname][randBlockEnum])
                    #/
                    
                    # classify gain/loss status for blocks
                    blocks_status_rand = classify_window_blocks_status(blocks_rand,sample,rname)
                    blocks_emb_status_rand = classify_window_blocks_status(blocks_rand_emb,'emb',rname)
                    #/
                    
                    # Classify gain/loss status for window
                    window_status_rand = classify_window_status(blocks_status_rand,blocks_emb_status_rand)
                    #/
                    
                    # Save
                    if window_status_rand == 'gain':
                        init(samples_rname_windows_gain_rand,sample,{})
                        init(samples_rname_windows_gain_rand[sample],rname,[])
                        samples_rname_windows_gain_rand[sample][rname].append(blocks)
                    elif window_status_rand == 'loss':
                        init(samples_rname_windows_loss_rand,sample,{})
                        init(samples_rname_windows_loss_rand[sample],rname,[])
                        samples_rname_windows_loss_rand[sample][rname].append(blocks)
                    #/
            ##/
            # Increase block window start
            blockStart += 1
            #/
#/

# Summarize genes affected in chained blocks
gtfs_to_scan = {'exon', 'pre_miRNA', 'snRNA', '5UTR', 'tRNA', 'ncRNA', 'start_codon', 'miRNA', 'snoRNA', 'rRNA', 'stop_codon'}
genes_sampleHits_gain = compile_sample_rname_block_dict(samples_rname_windows_gain,gtfs_to_scan=gtfs_to_scan)
genes_sampleHits_loss = compile_sample_rname_block_dict(samples_rname_windows_loss,gtfs_to_scan=gtfs_to_scan)
genes_sampleHits_gain_rand = compile_sample_rname_block_dict(samples_rname_windows_gain_rand,gtfs_to_scan=gtfs_to_scan)
genes_sampleHits_loss_rand = compile_sample_rname_block_dict(samples_rname_windows_loss_rand,gtfs_to_scan=gtfs_to_scan)
#/
# Transform format, from: "gene -> samples" to: "sample -> genes"
samples_genes_gain = {}
samples_genes_loss = {}
for gene,samples in genes_sampleHits_gain.items():
    for sample in samples:
        init(samples_genes_gain,sample,set())
        samples_genes_gain[sample].add(gene)
for gene,samples in genes_sampleHits_loss.items():
    for sample in samples:
        init(samples_genes_loss,sample,set())
        samples_genes_loss[sample].add(gene)
#/
# Check number of blocks found up/down for each sample
samples_numWindows_gain = {}
samples_numWindows_loss = {}
for sample in SAMPLES_noEmb:
    windows_gained = 0
    try:
        for rname,blocks in samples_rname_windows_gain[sample].items():
            windows_gained += len(blocks)
    except:
        pass
    
    windows_loss = 0
    try:
        for rname,blocks in samples_rname_windows_loss[sample].items():
            windows_loss += len(blocks)
    except:
        pass
    
    samples_numWindows_gain[sample] = windows_gained
    samples_numWindows_loss[sample] = windows_loss
      
if 0 and 'plot':
    plot_keyVal_arr(sorted(samples_numWindows_gain.items(),key=lambda x: x[0]),xticksrotate_deg=90,separateTimelines=True)
    plot_keyVal_arr(sorted(samples_numWindows_loss.items(),key=lambda x: x[0]),xticksrotate_deg=90,separateTimelines=True)
    
    if 0 and 'dump numbers to file':
        with open(output_dir+'/'+'windows_gained.tsv','w') as nf:
            writeArr = ['sample','gained_abus']
            nf.write('\t'.join(map(str,writeArr))+'\n')
            for sample,count in sorted(samples_numWindows_gain.items(),key=lambda x: x[0]):
                writeArr = [sample,count]
                nf.write('\t'.join(map(str,writeArr))+'\n')
        
        with open(output_dir+'/'+'windows_lost.tsv','w') as nf:
            writeArr = ['sample','lost_abus']
            nf.write('\t'.join(map(str,writeArr))+'\n')
            for sample,count in sorted(samples_numWindows_loss.items(),key=lambda x: x[0]):
                writeArr = [sample,count]
                nf.write('\t'.join(map(str,writeArr))+'\n')
#/
### INFO
if 0 and 'number of samples per gene':
    numSamples = []
    for gene,samples in genes_sampleHits_gain.items():
        numSamples.append(len(samples))
    plotHist(numSamples)
    
    numSamples_req = 8
    with open(output_dir+'/'+'genes_gain_winSize5_numSamples'+str(numSamples_req)+'plus.tsv','w') as nf:
        for gene,samples in genes_sampleHits_gain.items():
            if len(samples) >= numSamples_req:
                nf.write(gene+'\n')
                
if 0 and 'output genes gained/lost per sample':
    tmp_out = output_dir+'/'+'genes_gained_per_sample'
    mkdir(tmp_out)
    for sample,genes in samples_genes_gain.items():
        with open(tmp_out+'/'+sample+'.tsv','w') as nf:
            nf.write('\n'.join(genes))
            
    with open(tmp_out+'/ALL_gained.tsv','w') as nf:
        for sample,genes in samples_genes_gain.items():
            nf.write(sample+'\t'+'\t'.join(genes)+'\n')
            
    with open(tmp_out+'/ALL_lost.tsv','w') as nf:
        for sample,genes in samples_genes_loss.items():
            nf.write(sample+'\t'+'\t'.join(genes)+'\n')
            
if 0 and 'make boxplots for sample coverage ratio in all genes covered by >3 samples':
    # (to check if these genes tend to be up in all samples [artefact?] or is clearly up in few samples)
    genes_samples_vals = [] #boxplot y vals for samples
    genes_calledSamples_vals = {} # array per sample (marked with gain) of scatter y vals
    genes_samples_ageVals = {} # age-> gene_sample_vals (For scatterplot values for samples, grouped by age)
    genes_targetTimeline_vals = {} # sample -> vals_per_gene
    genes_samples_stdevs = []
    genes_emb_vals = [] #scatter plot y vals for emb
    genes_emb_stdevs = []
    xlabels = []
    for gene,samples in genes_sampleHits_gain.items():
        
        # Check if gene is part of selection (from DAVID enriched  clusters)
        if 1:
            skipGene = True
            # Stress response?, >3 samples
            geneScan = 'FBGN0013276, FBGN0013277, FBGN0001229, FBGN0001228, FBGN0001225, FBGN0001224, FBGN0001227, FBGN0266599, FBGN0001226, FBGN0001223, FBGN0001233, FBGN0001230, FBGN0051354'
            # Dorsal closure
            #geneScan = 'FBGN0001297, FBGN0283426, FBGN0243512, FBGN0002121, FBGN0000308, FBGN0000414, FBGN0000575, FBGN0004657'
            # unfolded protein binding
            #geneScan = 'FBGN0013277, FBGN0001223, FBGN0035982, FBGN0051354, FBGN0010621, FBGN0001229, FBGN0001225, FBGN0266599, FBGN0001227'
            # glutathione transferase activity
            #geneScan = 'FBGN0042206, FBGN0063496, FBGN0063495, FBGN0086348, FBGN0038020, FBGN0001149, FBGN0035904, FBGN0035907, FBGN0035906, FBGN0010038, FBGN0063494, FBGN0063493, FBGN0063492, FBGN0010039, FBGN0063491, FBGN0010226'
            # extracellular matrix structural constituent, "many samples", found by SNPs with AF>0.5
            #geneScan = 'FBGN0053300, FBGN0053196, FBGN0052580, FBGN0036181'
            # transport along microtubule, "one sample enriched", found by SNPs with AF>0.5
            #geneScan = 'FBGN0003654, FBGN0283433'

            for gene2 in geneScan.split(', '):
                if gene[4:] == gene2[4:]:
                   skipGene = False
            if skipGene: continue
        #/
        if len(samples) > 0:
            if len(genes_samples_vals) >= 20: break
            coverage_ratios = []
            stdevs = []
            age_ratios = {}
            
            gene_entry = GTF_idx_map[ gene_GTF_idx_map[gene] ]
            rname = gene_entry['rname']
            rcoords = gene_entry['rcoords']
            rstart_roundDownThousand = int(rcoords[0]/1000)*1000
            
            ## Parse cov ratios for all freec entries that covers gene
            samples_covRatios_pre = {}
            while True:
                freec_entry = rname_freecs[rname][rstart_roundDownThousand]
                numBP_ovlp = getRangeOvlp(freec_entry['rcoords'],gene_entry['rcoords'])
                # check if over full bin (or gene size)
                if numBP_ovlp >= min([gene_entry['rcoords'][-1]-gene_entry['rcoords'][0] , 1000-1]): # 999numbp is ovlp for 1000 binsize
                    for sample,cov_ratio in freec_entry['samples'].items():
                        init(samples_covRatios_pre,sample,[])
                        samples_covRatios_pre[sample].append(cov_ratio)
                #/
                # Check if 0 overlap, then break loop (we now slide outside gene region)
                if numBP_ovlp < 0: break
                #/
                rstart_roundDownThousand += 1000 #move to next bin
                
            samples_covRatios = {}
            for sample,cov_ratios in samples_covRatios_pre.items():
                samples_covRatios[sample] = mean(cov_ratios)
                
            if not samples_covRatios: continue #skip if we found no cov ratio
            ##/

            for sample in sorted(SAMPLES):
                cov_ratio = samples_covRatios[sample]
                
                # Save embryonal value
                if sample == 'emb':
                    genes_emb_vals.append(cov_ratio)
                    genes_emb_stdevs.append(sample_rname_stats[sample][rname]['stdev'])
                #/
                # Save other values
                else:
                    coverage_ratios.append(cov_ratio)
                    
                    for age in ('1','2','3',):
                        if SAMPLES_age[sample] == age:
                            init(age_ratios,age,[])
                            age_ratios[age].append(cov_ratio)
                            
                    
                    stdevs.append(sample_rname_stats[sample][rname]['stdev'])
                    
                    # Save values for samples used in gene copy number gain call to mark with scatterplot
                    init(genes_calledSamples_vals,sample,[])
                    if sample in samples:
                        genes_calledSamples_vals[sample].append(cov_ratio)
                    else:
                        genes_calledSamples_vals[sample].append(0)
                    #/
                #/
                # Check if we have target timeline, then save separately
                target_sample_timeline = 'D5-S-P5-203'
                if SAMPLES_timeline[target_sample_timeline] == SAMPLES_timeline[sample]:
                    init(genes_targetTimeline_vals,sample,[])
                    genes_targetTimeline_vals[sample].append(cov_ratio)
                #/
            
            # Save to outer
            genes_samples_vals.append(coverage_ratios)
            genes_samples_stdevs.append(stdevs)
            
            for age,vals in age_ratios.items():
                init(genes_samples_ageVals,age,[])
                genes_samples_ageVals[age].append(vals)
            
            xlabels.append(gene)
            #/
    
    fig = plt.figure(figsize=(20,11))
    plt.boxplot(genes_samples_vals,sym='')
    ax = plt.gca()
    
    # Embryo scatter
    scatter_xaxis = []
    for i in range(len(genes_samples_vals)):
        scatter_xaxis.append((i+1)-0.3)
    ax.scatter(scatter_xaxis,genes_emb_vals,facecolor='white',edgecolor='black',s=100)
    #/
    
    # Target timeline scatter
    if 1:
        for sample,vals in genes_targetTimeline_vals.items():
            tmp_xaxis = []
            for i in range(len(vals)):
                tmp_xaxis.append((i+1)+0.3)
                
            tmp_shape = None
            if SAMPLES_age[sample] == '1':  tmp_shape = 'v'
            if SAMPLES_age[sample] == '2':  tmp_shape = 'x'
            if SAMPLES_age[sample] == '3':  tmp_shape = '^'
            ax.scatter(tmp_xaxis,vals,color='black',marker=tmp_shape)
    #/
    
    # Age scatter
    scatter_agevals_xaxis = {} # age -> x_coordinate_offset per gene
    for age,genes_vals in genes_samples_ageVals.items():
        age_offset = 0
        if age == '1':   age_offset -= 0.2
        elif age == '3': age_offset += 0.2
        
        for i in range(len(genes_vals)):
            xaxis_val = (i+1)+age_offset
            init(scatter_agevals_xaxis,age,[])
            scatter_agevals_xaxis[age].append(xaxis_val)
        
    for age,genes_vals in genes_samples_ageVals.items():
        plot_col = 'None'
        if age == '1':
            plot_col = '#BCBEC0'
        elif age == '2':
            plot_col = '#6D6E71'
        elif age == '3':
            plot_col = '#000000'
        for xaxis_pos,gene_vals in enumerate(genes_vals):
            ax.scatter([scatter_agevals_xaxis[age][xaxis_pos]]*len(gene_vals),gene_vals,facecolor=plot_col,edgecolor='black',s=100)
    #/

    # Samples marked with copy number gain in gene
    if 0:
        for sample,vals in genes_calledSamples_vals.items():
            ax.scatter(scatter_xaxis,vals,color='red')
    #/

    ax.set_xticklabels(xlabels,rotation=90)
    ax.set_ylim(bottom=1)
    
    plt.tight_layout()
    
    if 0 and 'savefig?':
        plt.savefig(output_dir+'/'+'genes_scatterplot.pdf')
        
    plt.show()

###/
    
### Plot the distribution of window-gained-blocks for selected samples
# Grab values for bins that are in gained windows
numWindows_check_thresh = 100 #skip samples that do not have at least this number of windows
rnames_start_coords = {} # rname -> block_starts ==> keep track of start-coordinates added. will use later to fill zeros between these ranges
samples_rnames_gainWin_binVals = {} # fill as rname -> block_start -> status (1=gain?).
for sample,rname_windows in samples_rname_windows_gain.items():
    if not samples_numWindows_gain[sample] >= numWindows_check_thresh: continue #skip if too few blocks
    
    for rname,windows in rname_windows.items():
        for window in windows:
            # Add positions by blocks
            for block in window:
                block_start = block['rcoords'][0]
                
                init(samples_rnames_gainWin_binVals,sample,{})
                init(samples_rnames_gainWin_binVals[sample],rname,{})
                samples_rnames_gainWin_binVals[sample][rname][block_start] = 1
                
                init(rnames_start_coords,rname,[])
                rnames_start_coords[rname].append(block_start)
            #/
            
            # Add all positions in window
            if 1 and 'add all positions of a window (may have gaps due to repeatmasked block)':
                win_start = window[0]['rcoords'][0]
                win_end = window[-1]['rcoords'][0]
                
                for start_pos in range(win_start,win_end+1000,1000):
                    init(samples_rnames_gainWin_binVals,sample,{})
                    init(samples_rnames_gainWin_binVals[sample],rname,{})
                    
                    if not start_pos in samples_rnames_gainWin_binVals[sample][rname]:
                        samples_rnames_gainWin_binVals[sample][rname][start_pos] = 1
            #/
#/
# Fill zero's
if 1:
    for sample,rnames_block_starts in samples_rnames_gainWin_binVals.items():
        for rname,block_starts in rnames_block_starts.items():
            for i in range(0,max(rnames_start_coords[rname]),1000):
                if not i in samples_rnames_gainWin_binVals[sample][rname]:
                    samples_rnames_gainWin_binVals[sample][rname][i] = 0
#/
# Grab genes that we want to plot...
genes_to_plot = []
for i in 'FBGN0013276, FBGN0013277, FBGN0001229, FBGN0001228, FBGN0001225, FBGN0001224, FBGN0001227, FBGN0266599, FBGN0001226, FBGN0001223, FBGN0001233, FBGN0001230, FBGN0051354'.split(', '):
    genes_to_plot.append('FBgn'+i.lower().replace('fbgn',''))
    
rnames_genes_to_plot = {}
for gene in genes_to_plot:
    if gene in gene_GTF_idx_map:
        gene_entry = GTF_idx_map[gene_GTF_idx_map[gene]]
        init(rnames_genes_to_plot,gene_entry['rname'],[])
        rnames_genes_to_plot[gene_entry['rname']].append(gene_entry)
#/
# Plot distribution per sample and chromosome
if 0 and 'plot distribution per chromosome?':
    # Plot per sample (chromosomes are rows)
    chroms_enums = {'2L':0,'2R':1,'3L':2,'3R':3,'X':4}
    for sample,rnames_vals in sorted(samples_rnames_gainWin_binVals.items(), key=lambda x: x[0]):
        fig,ax = plt.subplots(5,1,figsize=(15,10))
        
        for rname,vals in sorted(rnames_vals.items(), key=lambda x: x[0]):
            if rname == '4': continue
            enum = chroms_enums[rname]
            
            x_vals = []
            y_vals = []
            
            for pos,val in sorted(vals.items(), key=lambda x: x[0]):
                x_vals.append(pos)
                y_vals.append(val)
            
            
            ax[enum].plot(x_vals,y_vals)
            ax[enum].set_ylabel(rname,rotation='vertical')
            #plt.title(sample + ' || ' + rname)
            #plt.xscale('log')
            
            # Check if plot genes
            if 1 and 'plot genes?':
                if rname in rnames_genes_to_plot:
                    for gene_entry in rnames_genes_to_plot[rname]:
                        ax[enum].scatter(gene_entry['rcoords'][0],1.1,color='red')
            #/
            
        ax[0].set_title(sample)
        plt.show()
    #/
    
    # Plot per chromosome (samples are rows)
    for rname_target in ('2L','2R','3L','3R','X',):
        fig,ax = plt.subplots(len(samples_rnames_gainWin_binVals),1,figsize=(15,10),sharex=True)
        
        for enum,(sample,rnames_vals) in enumerate(sorted(samples_rnames_gainWin_binVals.items(), key=lambda x: x[0])):
            for rname,vals in sorted(rnames_vals.items(), key=lambda x: x[0]):
                if not rname == rname_target: continue
            
                x_vals = []
                y_vals = []
                
                for pos,val in sorted(vals.items(), key=lambda x: x[0]):
                    x_vals.append(pos)
                    y_vals.append(val)
                    
                ax[enum].plot(x_vals,y_vals)
                ax[enum].set_ylabel(sample,rotation='vertical')
                
                # Check if plot genes
                if 1 and 'plot genes?':
                    if rname in rnames_genes_to_plot:
                        for gene_entry in rnames_genes_to_plot[rname]:
                            ax[enum].scatter(gene_entry['rcoords'][0],1.1,color='red')
                #/
                
                
        ax[0].set_title(rname_target)
        plt.show()
    #/
#/
    
# Get lengths of windows
samples_windows_lenDistri = {}
for sample,rnames_block_starts in samples_rnames_gainWin_binVals.items():
    for rname,block_starts in rnames_block_starts.items():
        # Chain block starts on current rname
        block_starts_chained = []
        for block_start,val in sorted(block_starts.items(),key=lambda x: x[0]):
            if val > 0: # ignore if value is 0
                # check if previous block is adjacent to current
                if block_starts_chained and (block_start - block_starts_chained[-1][-1] == 1000):
                    block_starts_chained[-1].append(block_start)
                
                # else, init new chain
                else:
                    block_starts_chained.append([block_start])
        #/
        # Save window lengths to outer
        init(samples_windows_lenDistri,sample,[])
        for chain in block_starts_chained:
            samples_windows_lenDistri[sample].append(len(chain))
        #/
#/
###/

### Get copy-number status of selected genes in other Dmel cell lines
# Parse copynumber file
lee_gene_CN_file = '/Disk/user/gene_copy_number_lee2014.tsv' # copy-number data sheet exported from Lee, et. al. (2014)
lee_gene_CNs = {}
with open(lee_gene_CN_file,'r') as f:
    header = None
    for enum,line in enumerate(f):
        line = line.strip('\n')
        line = line.split('\t')
        
        # parse header
        if enum == 0:
            header = line
            continue
        #/
        
        # parse gene data
        gene_id = line[0]
        if not gene_id in genes_to_plot: continue #skip if gene is not wanted
        
        lines_CN = {}
        for header_enum,header_entry in enumerate(header):
            if header_entry.find('CopyNumber') != -1:
                celline_ID = header_entry.replace('_CopyNumber','')
                celline_CN = line[header_enum]
                try:
                    celline_CN = int(celline_CN)
                except:
                    celline_CN = -1
                
                lines_CN[celline_ID] = int(celline_CN)
        #/
        # Save
        lee_gene_CNs[gene_id] = lines_CN
        #/
#/
# Enter cell line ploidies, from Lee "Ploidy of cell lines":
#   "We scored ... as minimally diploid ... tetraploid"
cell_line_ploidies = {}
cell_line_ploidies['1182-4H'] = 2
cell_line_ploidies['BG3-c2'] = 2
cell_line_ploidies['Cl.8-BG'] = 2
cell_line_ploidies['Cl.8-DM'] = 2
cell_line_ploidies['D16-c3'] = 4
cell_line_ploidies['D17-c3'] = 4
cell_line_ploidies['D20-c2'] = 2
cell_line_ploidies['D20-c5'] = 2
cell_line_ploidies['D4-c1'] = 2
cell_line_ploidies['D8'] = 2
cell_line_ploidies['D9'] = 2
cell_line_ploidies['Kc167'] = 4
cell_line_ploidies['L1'] = 2
cell_line_ploidies['mbn2'] = 4
cell_line_ploidies['S1'] = 2
cell_line_ploidies['S2-DRSC-DM'] = 4
cell_line_ploidies['S2-DRSC-BO'] = 4
cell_line_ploidies['S2R+'] = 4
cell_line_ploidies['S3'] = 4
cell_line_ploidies['Sg4'] = 4
cell_line_ploidies['W2'] = 2

#/
# Plot copy-number in selected genes
if 0 and 'plot?':
    fig,ax = plt.subplots(len(lee_gene_CNs),1,figsize=(15,15),sharex=True)
    x_labels = None
    x_lim = [-0.5,20.5]
    for enum,(gene_id,cellines_CNs) in enumerate(sorted(lee_gene_CNs.items(),key=lambda x: x[0])):
        x_vals = []
        y_vals = []
        for celline,CN in sorted(cellines_CNs.items(), key=lambda x: x[0]):
            CN_deviation = cell_line_ploidies[celline] - CN
            x_vals.append(celline)
            y_vals.append(CN_deviation)
            
        if not x_labels:        x_labels = x_vals
        
        ax[enum].plot(x_lim,[4,4],color='black',zorder=0,linewidth=0.5)
        ax[enum].plot(x_lim,[2,2],color='black',zorder=0,linewidth=0.5)
        ax[enum].plot(x_lim,[0,0],color='black',zorder=0,linewidth=1)
        ax[enum].plot(x_lim,[-2,-2],color='black',zorder=0,linewidth=0.5)
        ax[enum].plot(x_lim,[-4,-4],color='black',zorder=0,linewidth=0.5)
        ax[-1].set_xlim(x_lim)
        
        #ax[enum].grid(axis = 'y', zorder=0)
        ax[enum].bar(list(range(len(x_vals))),y_vals,zorder=3)
        ax[enum].set_ylabel(gene_id,rotation='vertical',fontsize=8)
        ax[enum].set_yticks([0,1,2,3,4,5,6])

        
        
    ax[0].set_title(rname_target)
    ax[-1].set_xticks(list(range(len(x_vals))))
    ax[-1].set_xticklabels(x_labels,rotation=90)
    plt.show()
#/
###/

### Check frequency of gained genes in lines vs. published data
# Parse data
files_dir = '/Disk/user/12864_2004_176_MOESM' # lists of stress-related genes from Girardot, et. al. (2004)
files_genes = {}
for file_ in os.listdir(files_dir):
    if file_.endswith('.txt'):
        with open(files_dir+'/'+file_,'r') as f:
            header_line = None
            for enum,line in enumerate(f):
                line = line.strip('\n')
                line = line.split('\t')
                # parse header line
                if enum == 0:
                    header_line = line
                    continue
                #/
                id_spot,gene_id,symbol,detect_pval,meanP15,meanP5,meanH1,meanT12 = line[:8]
                
                name,product,function = line[-3:]
                
                file_name = file_.split('.')[0]
                init(files_genes,file_name,{})
                files_genes[file_name][gene_id] = {'detect_pval':detect_pval, 'meanP15':meanP15,
                                                   'meanP5':meanP5,'meanH1':meanH1,
                                                   'meanT12':meanT12,
                                                   'name':name,'product':product,
                                                   'function':function}
files_genes_pan = set()
for _,genes in files_genes.items():
    files_genes_pan.update(genes)
#/
# Calc the expected percentage of hits to each category (fileName)
files_genes_rate = {}
for fileName,genes in files_genes.items():
    numGenes = len(genes)
    percGenes = round(len(genes) / len(gene_GTF_idx_map),4)
    internal_percGenes = round(len(genes)/len(files_genes_pan),4)
    
    files_genes_rate[fileName] = {'numGenes':numGenes,'percGenes':percGenes,
                    'internal_percGenes':internal_percGenes}
#/
# Get genes that are gained in samples and that are found in published data
samples_gainGenes_matched = {}
for sample,sample_genes in samples_genes_gain.items():
    if len(sample_genes) < 50: continue #skip samples with few hits to genes
    
    for file_name,file_genes in files_genes.items():
        # init file entry at sample
        init(samples_gainGenes_matched,sample,{})
        init(samples_gainGenes_matched[sample],file_name,{})
        #/
        for file_gene in file_genes:
            if file_gene in sample_genes:
                samples_gainGenes_matched[sample][file_name][file_gene] = file_genes[file_gene]
#/
# Check percentage of matched genes and relate it to the expected
samples_gainGenes_matched_rate = {}
for sample,fileNames_genes in samples_gainGenes_matched.items():
    for fileName,genes in fileNames_genes.items():
        numGenes = len(genes)
        percGenes = round(len(genes) / len(samples_genes_gain[sample]),4)
        init(samples_gainGenes_matched_rate,sample,{})
        samples_gainGenes_matched_rate[sample][fileName] = {'numGenes':numGenes,'percGenes':percGenes}
#/
# Dump for plot
if 0 and 'dump for plot?':
    tmp_out = output_dir+'/'+'genes_stress_matched_rate'
    mkdir(tmp_out)
    with open(tmp_out+'/'+'matched_rate_nums.tsv','w') as nf:
        for sample in samples_gainGenes_matched_rate:
            for fileName in samples_gainGenes_matched_rate[sample]:
                entry = 'numGenes'
                value = samples_gainGenes_matched_rate[sample][fileName][entry]
                writeArr = [sample,fileName,entry,value]
                nf.write('\t'.join(map(str,writeArr))+'\n')
#/
###/
