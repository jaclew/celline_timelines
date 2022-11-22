
import sys
import os
import matplotlib.pyplot as plt
from statistics import mean,median,stdev

from functions import * # imports from functions.py (must be in the same directory as current script)

### INPUTS
dir_of_dirs = '/Disk/user/jitterbug/runs' #path to where output directories of jitterbug is
jitterbug_basename = '<sample>.TE_insertions_paired_clusters.gff3' # will replace <sample> with folder names (of samples) in "dir_of_dirs"
###/INPUTS

### OUTPUTS
output_dir = '/Disk/user/jitterbug_analysis_output'
###/

### GLOBALS
global_vars_modules.global_idx_iterator = 0 # keeps track of added index via assignIdx function
chroms = ['2L','2R','3L','3R','X','4']
autosomes = ['2L','2R','3L','3R']
###/GLOBALS

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
    
### Parse GFF files
TE_idx_map = {}
SAMPLES = {}
for file_ in os.listdir(dir_of_dirs):
    sample = file_
    path1 = dir_of_dirs+'/'+sample+'/'+jitterbug_basename.replace('<sample>',sample)
    if not os.path.exists(path1): continue #skip if we did not find jitterbug output file in current directory
    
    with open(path1,'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split('\t')
            rname = line[0]
            rstart = int(line[3])
            rend = int(line[4])
            if rend == rstart:      rend+=1 # adjust rend position if same as rstart
            
            data = line[-1].split('; ')
            supp1 = int(data[0].split('=')[1])
            supp2 = int(data[1].split('=')[1])
            cluster_pair_id = int(data[2].split('=')[1])
            TEs1 = data[4].split('=')[1].split(', ')
            TEs2 = data[5].split('=')[1].split(', ')
            
            tmp_save = {'rname':rname,'rcoords':[rstart,rend],'supp1':supp1,'supp2':supp2,'TEs1':TEs1,'TEs2':TEs2,
                        'sample':sample,'cluster_pair_id':cluster_pair_id}
            assignIdx(tmp_save)
            TE_idx_map[tmp_save['idx']] = tmp_save
            
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
samples_enums = {}
for enum,sample in enumerate(list(sorted(SAMPLES_noEmb))):
    samples_enums[sample] = enum
###/

### Ban all TE calls which ovlp emb calls (using ovlpRanges function to find overlaps)
# setup ranges
fatten_margin = 1000
rname_rRanges = {}
for TE_idx,TE in TE_idx_map.items():
    rname = TE['rname']
    rcoords_wMarg = []
    init(rname_rRanges,rname,[])
    rname_rRanges[rname].append([TE['rcoords'][0]-fatten_margin,TE_idx,TE['rcoords'][-1]+fatten_margin])
#/
# find ovlps and find which TEs to ban
banned_TEs = set()
for rname,rangees in rname_rRanges.items():
    ovlps = computeOvlpRanges_wIdxs(rangees)
    for oS,TE_idxs,oE in ovlps:
        # check if embryo exist at position
        embExist = False
        for TE_idx in TE_idxs:
            TE = TE_idx_map[TE_idx]
            if TE['sample'] == 'emb':
                embExist = True
                break
        #/
        # If emb existed, then ban all other TE idxs at ovlpRange
        if embExist:
            for TE_idx in TE_idxs:
                banned_TEs.add(TE_idx)
        #/
#/
# Add marker at TEs if banned
INFO_embOvlp = {'tot':0,'yes':0,'isEmb':0}
for TE_idx,TE in TE_idx_map.items():
    TE['embOvlp'] = False
    INFO_embOvlp['tot'] += 1
    if TE['sample'] == 'emb':       INFO_embOvlp['isEmb'] += 1
    if TE_idx in banned_TEs:
        TE['embOvlp'] = True
        INFO_embOvlp['yes'] += 1
#/
###/

### INFO
if 0 and 'INFO plot abundance of TEs in samples':
    samples_vals = {}
    for TE in TE_idx_map.values():
        if TE['embOvlp'] and not TE['sample'] == 'emb': continue
        sample = TE['sample']
        init(samples_vals,sample,0)
        samples_vals[sample] += 1
        
    #if 'emb' in samples_vals:           del samples_vals['emb']
    vals_nSorted = sorted(samples_vals.items(),key=lambda x: x[0])
    vals_ageSorted = sorted(samples_vals.items(),key=lambda x: int(x[0].split('-')[-1]))
    plot_keyVal_arr(vals_nSorted,title='number of TEs',xticksrotate_deg=90,separateTimelines=True)
###/
### Count TE families (1) overall and (2) per sample
# numSupp for all, and count of all TEs
numSupps = []
TE_counts = {}
for TE_idx,TE in TE_idx_map.items():
    if not TE['rname'] in chroms: continue
    if TE['embOvlp']: continue

    numSupps.append(TE['supp1'])
    numSupps.append(TE['supp2'])
    
    for te in TE['TEs1']+TE['TEs2']:
        te = te.split('{}')[0]
        init(TE_counts,te,0)
        TE_counts[te] += 1

if 0 and 'plot?':
    plotHist(numSupps)

TE_counts_filt = []
for TE,count in sorted(TE_counts.items(),key=lambda x: x[1],reverse=True):
    if count < 100: continue
    TE_counts_filt.append([TE,count])
#/

# Counts of filtered TEs, per sample
sample_TE_counts = {}
for TE_idx,TE in TE_idx_map.items():
    if not TE['rname'] in chroms: continue
    if TE['embOvlp']: continue

    for te in TE['TEs1']+TE['TEs2']:
        te = te.split('{}')[0]
        if te in dict(TE_counts_filt):
            sample = TE['sample']
            init(sample_TE_counts,sample,{})
            init(sample_TE_counts[sample],te,0)
            sample_TE_counts[sample][te] += 1
#/
# Filter above by abundance
sample_TE_counts_filt = {}
for sample,tes_counts in sample_TE_counts.items():
    for te,count in tes_counts.items():
        if count < 50: continue
        init(sample_TE_counts_filt,sample,{})
        sample_TE_counts_filt[sample][te] = count
#/
# Find top X elements per sample
sample_TE_counts_topX_perc = {}
for sample,tes_counts in sample_TE_counts.items():
    for enum,(te,count) in enumerate(sorted(tes_counts.items(),key=lambda x: x[1],reverse=True)):
        if enum == 5: break
        te_abu_frac = round(count / sum(tes_counts.values()),4)
        init(sample_TE_counts_topX_perc,sample,{})
        sample_TE_counts_topX_perc[sample][te] = te_abu_frac
#/
###/

### Find which samples exist at ovlpRange regions
# setup ranges
fatten_margin = 1000
rname_rRanges = {}
for TE_idx,TE in TE_idx_map.items():
    rname = TE['rname']
    if not rname in chroms: continue
    rcoords_wMarg = []
    init(rname_rRanges,rname,[])
    rname_rRanges[rname].append([TE['rcoords'][0]-fatten_margin,TE_idx,TE['rcoords'][-1]+fatten_margin])
#/
# find ovlps and if an embro-overlapping call is made, then skip
TE_rname_ovlps = {}
for rname,rangees in rname_rRanges.items():
    ovlps = computeOvlpRanges_wIdxs(rangees)
    for oS,TE_idxs,oE in ovlps:
        # check if banned call (e.g., overlaps embryo) exist at position
        bannedExist = False
        for TE_idx in TE_idxs:
            TE = TE_idx_map[TE_idx]
            if TE['embOvlp']:
                bannedExist = True
                break
        if bannedExist: continue
        #/
        
        # Save ovlp
        init(TE_rname_ovlps,rname,[])
        if (not TE_rname_ovlps[rname]) or TE_rname_ovlps[rname][-1][-1] != oS: # check if init new
            TE_rname_ovlps[rname].append([oS,{'rname':rname,'TE_idxs':set(TE_idxs)},oE])
        else: #else update/extend previous entry
            TE_rname_ovlps[rname][-1][1]['TE_idxs'].update(TE_idxs)
            TE_rname_ovlps[rname][-1][-1] = oE
        #/
#/

# At multi-sample ovlpRanges, check what types of samples exist (i.e. is it intra-timeline? intra-age?)
timelines_allPresent_onlyTimeline = {} #timeline -> count
timelines_allPresent = {}
timelines_morePresent = {}
timelines_morePresent_onlyTimeline = {}
shared_loci_assignments_abu = {}

timelineUniq_anyTimelineSamples = {} # TE's in only one timeline, may appear in multiple samples intra-timeline
for rname,ovlpRangees in TE_rname_ovlps.items():
    for ovlpRangee in ovlpRangees:
        samples = set()
        for TE_idx in ovlpRangee[1]['TE_idxs']:
            samples.add(TE_idx_map[TE_idx]['sample'])
            
        #if len(samples) < 4: continue
        #if len(samples) > 15: continue
            
        # Count number of timeline samples
        rangee_timelines = {}
        for sample in samples:
            timeline = SAMPLES_timeline[sample]
            init(rangee_timelines,timeline,set())
            rangee_timelines[timeline].add(sample)
        
        if 0 and 'purge low-frequent timelines?':
            purge_timelines = []
            for timeline,samplees in rangee_timelines.items():
                if len(samplees) <= 1:
                    purge_timelines.append(timeline)
            for purge_timeline in purge_timelines:
                del rangee_timelines[purge_timeline]
        #/
        
        ## Save TE count where only one sample (or samples within that timeline) has the call
        if len(rangee_timelines) == 1:
            for sample in samples:
                init(timelineUniq_anyTimelineSamples,sample,[])
                
                # check so we didnt add ovlpRange chain previously
                saveEntry = True
                for existEntry in timelineUniq_anyTimelineSamples[sample]:
                    if existEntry['rname'] == rname and getRangeOvlp(existEntry['rcoords'],[ovlpRangee[0],ovlpRangee[-1]]) >= 0:
                        saveEntry = False
                        print(sample,'hey')
                #/
                
                if saveEntry:
                    for TE_idx in ovlpRangee[1]['TE_idxs']:
                        if TE_idx_map[TE_idx]['sample'] == sample:
                            timelineUniq_anyTimelineSamples[sample].append({'rname':rname,'rcoords':[ovlpRangee[0],ovlpRangee[-1]],'TE_idx':TE_idx})
        ##/
        
        # Count number of age samples
        rangee_ages = {}
        for sample in samples:
            age = SAMPLES_age[sample]
            init(rangee_ages,age,set())
            rangee_ages[age].add(sample)
        
        purge_ages = []
        for age,samplees in rangee_ages.items():
            if len(samplees) <= 1:
                purge_ages.append(age)
        for purge_age in purge_ages:
            del rangee_ages[purge_age]
        #/
        
        # Check which samples got what assignment
        samples_assignments = {}
        for timeline,samplees in rangee_timelines.items():
            if len(samplees) >= 2:
                for sample in samplees:
                    init(samples_assignments,sample,set())
                    samples_assignments[sample].add('timeline')
        for age,samplees in rangee_ages.items():
            if len(samplees) >= 2:
                for sample in samplees:
                    init(samples_assignments,sample,set())
                    samples_assignments[sample].add('age'+age)
        
        for sample in samples:
            if not sample in samples_assignments:
                samples_assignments[sample] = 'orphan'
        #/
        
        # Check if any timeline has all samples present
        for timeline,samplees in rangee_timelines.items():
            if samplees.intersection(TIMELINES_samples[timeline]) == TIMELINES_samples[timeline]:
                # check if only this timeline
                if len(TIMELINES_samples[timeline]) == len(samples):
                    init(timelines_allPresent_onlyTimeline,timeline,0)
                    timelines_allPresent_onlyTimeline[timeline] += 1
                
                # save
                init(timelines_allPresent,timeline,0)
                timelines_allPresent[timeline] += 1
            
            # save if multiple samples in timeline were present
            init(timelines_morePresent,timeline,0)
            timelines_morePresent[timeline] += 1
            
            # save if only these timeline samples were present
            if len(samples) == len(samplees):
                init(timelines_morePresent_onlyTimeline,timeline,0)
                timelines_morePresent_onlyTimeline[timeline] += 1
        #/
        
        # Compile assignments
        assignments_samples = {}
        for sample,assignments in samples_assignments.items():
            for assignment in assignments:
                if assignment.find('age') != -1:
                    init(assignments_samples,assignment,set())
                    assignments_samples[assignment].add(sample)
                elif assignment == 'timeline':
                    timeline = SAMPLES_timeline[sample]
                    init(assignments_samples,timeline,set())
                    assignments_samples[timeline].add(sample)
        #/
        
        for assignment,samples in assignments_samples.items():
            if assignment.find('age') == -1:    assignment += '_num='+str(len(samples))
            if assignment.find('age') != -1:
                if len(samples) < 4: continue
            init(shared_loci_assignments_abu,assignment,0)
            shared_loci_assignments_abu[assignment] += 1
#/


## Transform entries in "timelineUniq_anyTimelineSamples" to counts
if 0 and 'plot?':
    timelineUniq_anyTimelineSamples_counts = {}
    for sample,entries in timelineUniq_anyTimelineSamples.items():
        timelineUniq_anyTimelineSamples_counts[sample] = len(entries)
    
    vals_nSorted = sorted(timelineUniq_anyTimelineSamples_counts.items(),key=lambda x: x[0])
    vals_ageSorted = sorted(timelineUniq_anyTimelineSamples_counts.items(),key=lambda x: int(x[0].split('-')[-1]))
    plot_keyVal_arr(vals_nSorted,title='number of TEs',xticksrotate_deg=90,separateTimelines=True)
    
    ## Get most abundant TEs across all samples and per sample
    TE_types_abus = {}
    sample_TE_types_abus = {}
    for sample,entries in timelineUniq_anyTimelineSamples.items():
       for entry in entries:
           TE = TE_idx_map[entry['TE_idx']]
           for strint in ('1','2',):
               TE_IDs = TE['TEs'+strint]
               TE_IDs_filt = set()
               for TE_ID in TE_IDs:
                   TE_IDs_filt.add(TE_ID.split('{}')[0])

               for TE_ID in TE_IDs_filt:
                   init(TE_types_abus,TE_ID,0)
                   TE_types_abus[TE_ID] += 1

                   init(sample_TE_types_abus,sample,{})
                   init(sample_TE_types_abus[sample],TE_ID,0)
                   sample_TE_types_abus[sample][TE_ID] += 1

    TE_types_abus_top20 = {}
    for TE_ID,abu in sorted(TE_types_abus.items(),key=lambda x: x[1],reverse=True):
        TE_types_abus_top20[TE_ID] = abu
        if len(TE_types_abus_top20) >= 20:
            break
    ##/

    
    if 0 and 'dump numbers to file':
        with open(output_dir+'/'+'TEs_abus.tsv','w') as nf:
            writeArr = ['sample','TE count']
            nf.write('\t'.join(map(str,writeArr))+'\n')
            for sample,count in sorted(timelineUniq_anyTimelineSamples_counts.items(),key=lambda x: x[0]):
                writeArr = [sample,count]
                nf.write('\t'.join(map(str,writeArr))+'\n')
                
##/
###/
