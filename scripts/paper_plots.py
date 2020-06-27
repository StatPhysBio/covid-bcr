import string
import sys
from utils import *
import matplotlib
from matplotlib import rcParams
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

colors = {"Severe": "#D55E00",
          'Mild': "#009E73",
          'Moderate':"#E69F00",
          'Healthy':"#0072B2"}
lightercolors = {"Severe": (255/255, 138/255, 47/255),
          'Mild': (46/255, 255/255, 200/255),
          'Moderate': (255/255, 212/255, 110/255),
          'Healthy': (144/255, 215/255, 255/255)}
observables = ['v gene', 'j gene', 'cdr3 length', 'vd ins', 'dj ins', 'vd del', 'dj del']
inbinning = {'cdr3 length': np.linspace(2,36,18),'vd ins':np.linspace(0,30,31), 'dj ins':np.linspace(0,30,31),
           'vd del':np.linspace(0,35,36), 'dj del':np.linspace(0,42,43)}

def geo_mean(iterable):
    return np.exp(np.mean(np.log(iterable)))

def geo_std(iterable):
    return np.exp(np.std(np.log(iterable)))

def cohort_values(counters,observable,binning,indict=False, geo=False):
    if not geo:
        avg_func = np.mean
        std_func = np.std
    else:
        avg_func = geo_mean
        std_func = geo_std
    if indict:
        std = {}
        mean = {}        
        totals = [sum(c.values()) for c in counters]
        x = []
        for key in counters[0]:
            mean[key] = avg_func([c[key]/totals[i] for i, c in enumerate(counters)])
            std[key] = std_func([c[key]/totals[i] for i, c in enumerate(counters)])            
        avgs = np.array(list(mean.values()))
        stds = np.array(list(std.values()))
        x = list(mean.keys())
    else:
        x = binning[observable][:-1]
        std, mean, freqs = [], [], []
        for c in counters:
            freq, bins = np.histogram(c, binning[observable],density=True)
            freqs.append(freq)
        for i,value in enumerate(freqs[0]):
            mean.append(avg_func([c[i] for c in freqs]))
            std.append(std_func([c[i] for c in freqs]))
        avgs = np.array(mean)
        stds = np.array(std)
    where_are_NaNs = np.isnan(avgs)
    avgs[where_are_NaNs] = 0
    where_are_NaNs = np.isnan(stds)
    stds[where_are_NaNs] = 0
    return avgs, stds, x

def cohort_plot(cohort_dict, key,in_ax=None,geo=False):
    if in_ax is None:
        fig, ax = plt.subplots(figsize=(4, 4),dpi=300)
    else:
        ax = in_ax
    for s in cohort_dict:
        x = cohort_dict[s][key][-1]
        ax.plot(x,
                cohort_dict[s][key][0],
                alpha=1.0,
                color=colors[s], label=s, linewidth = 1.0)
        if geo:
            ax.fill_between(x, 
                            cohort_dict[s][key][0]/cohort_dict[s][key][1], 
                            cohort_dict[s][key][0]*cohort_dict[s][key][1], 
                            color=colors[s], alpha=0.15, linewidth=0)
        else:
            ax.fill_between(x, 
                            cohort_dict[s][key][0]-cohort_dict[s][key][1], 
                            cohort_dict[s][key][0]+cohort_dict[s][key][1], 
                            color=colors[s], alpha=0.15, linewidth=0)
    if in_ax is None:
        return fig,ax
    else:
        return ax
    
def format_axes(in_ax, labelsize=15, ticksize=12, legendsize=7):
    in_ax.xaxis.set_minor_locator(AutoMinorLocator())
    in_ax.yaxis.set_minor_locator(AutoMinorLocator())
    in_ax.xaxis.label.set_size(labelsize)
    in_ax.yaxis.label.set_size(labelsize)
    in_ax.tick_params(axis='both', which='both', direction='in', width=1,
                   bottom=True, top=True, left=True, right=True)
    in_ax.tick_params(axis='both', which='major', direction='in',length=4,labelsize=ticksize,color='grey',
                   bottom=True, top=True, left=True, right=True)
    in_ax.tick_params(axis='both', which='minor', direction='in',length=2, color='grey',
                   bottom=True, top=True, left=True, right=True)
    in_ax.legend(fontsize=legendsize)

def make_cdr3_plot(to_plot, in_ax,geo=False, labelsize=15, ticksize=12, legendsize=7):
    cohort_plot(to_plot, 'cdr3 length',in_ax=in_ax,geo=geo)
    in_ax.set_xlabel('CDR3 length [aa]',fontname="Arial")
    in_ax.set_ylabel("PDF",fontname="Arial")
    format_axes(in_ax, labelsize=labelsize, ticksize=ticksize, legendsize=legendsize)
    
def make_vddel_plot(to_plot, in_ax,geo=False, labelsize=15, ticksize=12, legendsize=7):
    cohort_plot(to_plot, 'vd del',in_ax=in_ax,geo=geo)
    in_ax.set_xlabel('VD deletions [nt]',fontname="Arial")
    in_ax.set_ylabel("PDF",fontname="Arial")
    format_axes(in_ax,labelsize=labelsize, ticksize=ticksize, legendsize=legendsize)
    
def make_vdins_plot(to_plot, in_ax,geo=False, labelsize=15, ticksize=12, legendsize=7):
    cohort_plot(to_plot, 'vd ins',in_ax=in_ax,geo=geo)
    in_ax.set_xlabel('VD insertions [nt]',fontname="Arial")
    in_ax.set_ylabel("PDF",fontname="Arial")
    format_axes(in_ax, labelsize=labelsize, ticksize=ticksize, legendsize=legendsize)
    
def make_djdel_plot(to_plot, in_ax,geo=False, labelsize=15, ticksize=12, legendsize=7):
    cohort_plot(to_plot, 'dj del',in_ax=in_ax,geo=geo)
    in_ax.set_xlabel('DJ deletions [nt]',fontname="Arial")
    in_ax.set_ylabel("PDF",fontname="Arial")
    format_axes(in_ax, labelsize=labelsize, ticksize=ticksize, legendsize=legendsize)
    
def make_djins_plot(to_plot, in_ax,geo=False, labelsize=15, ticksize=12, legendsize=7):
    cohort_plot(to_plot, 'dj ins',in_ax=in_ax,geo=geo)
    in_ax.set_xlabel('DJ insertions [nt]',fontname="Arial")
    in_ax.set_ylabel("PDF",fontname="Arial")
    format_axes(in_ax, labelsize=labelsize, ticksize=ticksize, legendsize=legendsize)

def cohort_bar(plot_dict, in_cohort_dict, key, yaxis_upper=None,
               ax=None, labelsize=10, ticksize=8, legendsize=6):
    if ax is None:
        fig = plt.figure(dpi=300,figsize=(16,4))
        ax = fig.add_subplot(111)
    width = 1.0 / len(plot_dict) - 0.01
    bars = []
    sorted_genes = [gene for _,gene in sorted(zip(plot_dict['Healthy'][key][0], plot_dict['Healthy'][key][-1]),
                                              key=lambda pair: pair[0], reverse=True)]

    if 'v' in key:
        good_genes = []
        for g in sorted_genes:
            bools = 0
            for s in plot_dict:
                idx = plot_dict[s][key][-1].index(g)
                value = plot_dict[s][key][0][idx]
                bools += value < 0.01
            if bools != 4:
                good_genes.append(g)
    else:
        good_genes = sorted_genes

    xlabels = good_genes
    default_x = np.arange(0,len(good_genes),1)
    xs,ys = [], []
        
    for i,s in enumerate(plot_dict):
        x = default_x + width*i 
        xs.append(x)
        indices = [plot_dict[s][key][-1].index(g) for g in good_genes]
        y = [plot_dict[s][key][0][k] for k in indices]
        ys.append(y)

    if key == 'v gene':
        ax.set_xlim(-0.5,xs[3][-1]+0.5)
        if yaxis_upper is not None:
            ax.set_ylim(0,yaxis_upper)
        else:
            ax.set_ylim(0,0.3)
    else:
        ax.set_xlim(-0.3,xs[3][-1]+0.3)
        if yaxis_upper is not None:
            ax.set_ylim(0,yaxis_upper)
        else:
            ax.set_ylim(0,0.55)
    middle = np.mean(np.arange(0,4)) - 0.5
    middle_xticks = default_x + middle*width
    ax.set_xticks(middle_xticks)
    ax.set_xticklabels(xlabels, fontsize=10, rotation=90,fontname="Arial")
    
        
    for i,s in enumerate(plot_dict):
        bar = ax.bar(xs[i],ys[i],width,#yerr=errs[i],
                      color=colors[s],label=s)
        for d in in_cohort_dict[s][key]:
            for gidx,rect in enumerate(bar):
                ax.scatter(xs[i][gidx],d[good_genes[gidx]],marker='.',
                           color=lightercolors[s],zorder=3,s=15,
                           edgecolors='black',linewidths=0.3)
    ax.set_ylabel('PDF',fontsize=15,fontname="Arial")
    
    format_axes(ax, labelsize=labelsize, ticksize=ticksize, legendsize=legendsize)
    

def get_seq_stats(in_files, geo=False):
    stats = {}
    for f in in_files:
        instats = unpickle(f)
        pat = list(instats.keys())[0]
        stats[pat] = {'v gene':{}, 'j gene':{}, 'cdr3 length': {},
                                          'vd ins':{}, 'dj ins': {},
                                          'vd del':{}, 'dj del': {}}
        for key in stats[pat]:
            for item in instats[pat]['in frame no indels'][key]:
                stats[pat][key][item] = len(instats[pat]['in frame no indels'][key][item])
                
    severity = {"Healthy":[],"Mild":[],"Moderate":[],"Severe":[]}
    cohort_data = {"Healthy":{},"Mild":{},"Moderate":{},"Severe":{}}
    to_plot = {"Healthy":{},"Mild":{},"Moderate":{},"Severe":{}}
    for key in CONST_DATA_DICT:
        s = CONST_DATA_DICT[key]['severity']
        for o in observables:
            cohort_data[s][o] = []
            to_plot[s][o] = []
        severity[s].append(stats[key])

    for s in severity:
        for key in severity[s][0]:
            counters = [count[key] for count in severity[s]]
            equalize_counters(counters, 1)
            for c in counters:
                if key == 'v gene' or key == 'j gene':
                    cohort_data[s][key].append(sort_dict_by_value(c))
                else:
                    datalist = []
                    for value in c:
                        datalist += [value]*c[value]
                    if key == 'cdr3 length':
                        datalist = np.array(datalist) / 3
                    cohort_data[s][key].append(datalist)
            if key == 'v gene' or key == 'j gene':
                continue
            else:
                to_plot[s][key] = cohort_values(cohort_data[s][key],key,inbinning,indict=False,geo=False)
            
    equalize_counters([d for s in cohort_data for d in cohort_data[s]['v gene']],1)
    equalize_counters([d for s in cohort_data for d in cohort_data[s]['j gene']],1)
    for s in severity:
        for key in ['v gene', 'j gene']:
            for c in cohort_data[s][key]:
                normalize_counter(c)
    for s in severity:
        for key in severity[s][0]:
            if key == 'v gene' or key == 'j gene':
                to_plot[s][key] = cohort_values(cohort_data[s][key],key,inbinning,indict=True,geo=False)
    return cohort_data, to_plot

def get_lin_stats(in_files, geo=False):
    stats = {}
    for f in in_files:
        patient = f.split("/")[-1].split("_")[0]
        instats = unpickle(f)
        for key in instats:
            instats[key] = sort_dict_by_value(Counter(instats[key]))
        del instats['d gene']
        stats[patient] = instats


    severity = {"Healthy":[],"Mild":[],"Moderate":[],"Severe":[]}
    cohort_data = {"Healthy":{},"Mild":{},"Moderate":{},"Severe":{}}
    to_plot = {"Healthy":{},"Mild":{},"Moderate":{},"Severe":{}}
    for key in CONST_DATA_DICT:
        s = CONST_DATA_DICT[key]['severity']
        for o in observables:
            cohort_data[s][o] = []
            to_plot[s][o] = []
        severity[s].append(stats[key])

    for s in severity:
        for key in severity[s][0]:
            counters = [count[key] for count in severity[s]]
            equalize_counters(counters, 1)
            for c in counters:
                if key == 'v gene' or key == 'j gene':
                    cohort_data[s][key].append(sort_dict_by_value(c))
                else:
                    datalist = []
                    for value in c:
                        datalist += [value]*c[value]
                    if key == 'cdr3 length':
                        datalist = np.array(datalist) / 3
                    cohort_data[s][key].append(datalist)
            if key == 'v gene' or key == 'j gene':
                continue
            else:
                to_plot[s][key] = cohort_values(cohort_data[s][key],key,inbinning,indict=False,geo=False)

    equalize_counters([d for s in cohort_data for d in cohort_data[s]['v gene']],1)
    equalize_counters([d for s in cohort_data for d in cohort_data[s]['j gene']],1)
    for s in severity:
        for key in ['v gene', 'j gene']:
            for c in cohort_data[s][key]:
                normalize_counter(c)
    for s in severity:
        for key in severity[s][0]:
            if key == 'v gene' or key == 'j gene':
                to_plot[s][key] = cohort_values(cohort_data[s][key],key,inbinning,indict=True,geo=False)
    return cohort_data, to_plot

lin_files = get_files('/gscratch/stf/zachmon/covid/stats','_conservative_lineage_stats_filtered.pickle')
lin_cd, lin_tp = get_lin_stats(lin_files)

fig8 = plt.figure(dpi=600,figsize=(10,5),constrained_layout=True)
gs1 = fig8.add_gridspec(nrows=2, ncols=3, wspace=0.025, hspace=0.05)
ax0 = fig8.add_subplot(gs1[:-1, :])
ax1 = fig8.add_subplot(gs1[-1, 0])
ax2 = fig8.add_subplot(gs1[-1, 1])
ax3 = fig8.add_subplot(gs1[-1, -1])
cohort_bar(lin_tp, lin_cd,'v gene',ax=ax0, yaxis_upper=0.21, labelsize=12, ticksize=10, legendsize=8)
make_cdr3_plot(lin_tp, ax1, labelsize=12, ticksize=10, legendsize=8)
make_vddel_plot(lin_tp, ax2, labelsize=12, ticksize=10, legendsize=8)
make_vdins_plot(lin_tp, ax3, labelsize=12, ticksize=10, legendsize=8)
for n, ax in enumerate([ax0, ax1, ax2, ax3]):
    ax.text(-0.1, 1.1, "("+string.ascii_uppercase[n]+")", transform=ax.transAxes, 
            size=12, weight='bold')
fig8.savefig('/gscratch/stf/zachmon/covid/lineages_greaterthan2_arithmetic_avg.pdf')

seq_files = get_files('/gscratch/stf/zachmon/covid/stats/','_stats_no_sings.pickle')
seq_cd, seq_tp = get_seq_stats(seq_files)

fig8 = plt.figure(dpi=600,figsize=(10,5),constrained_layout=True)
gs1 = fig8.add_gridspec(nrows=2, ncols=3, wspace=0.025, hspace=0.05)
ax0 = fig8.add_subplot(gs1[:-1, :])
ax1 = fig8.add_subplot(gs1[-1, 0])
ax2 = fig8.add_subplot(gs1[-1, 1])
ax3 = fig8.add_subplot(gs1[-1, -1])
cohort_bar(seq_tp, seq_cd,'v gene',ax=ax0, yaxis_upper=0.3, labelsize=12, ticksize=10, legendsize=8)
make_cdr3_plot(seq_tp, ax1, labelsize=12, ticksize=10, legendsize=8)
make_vddel_plot(seq_tp, ax2, labelsize=12, ticksize=10, legendsize=8)
make_vdins_plot(seq_tp, ax3, labelsize=12, ticksize=10, legendsize=8)
for n, ax in enumerate([ax0, ax1, ax2, ax3]):
    ax.text(-0.1, 1.1, "("+string.ascii_uppercase[n]+")", transform=ax.transAxes, 
            size=12, weight='bold')
fig8.savefig('/gscratch/stf/zachmon/covid/sequences_nosingletons_arithmetic_avg.pdf')

fig = plt.figure(dpi=600,figsize=(10,5), constrained_layout=True)
gs = fig.add_gridspec(nrows=2, ncols=4, wspace=0.025, hspace=0.05)
jax0 = fig.add_subplot(gs[0, 0:2])
jax1 = fig.add_subplot(gs[0, 2:4])
jax2 = fig.add_subplot(gs[1, 0])
jax3 = fig.add_subplot(gs[1, 1])
jax4 = fig.add_subplot(gs[1, 2])
jax5 = fig.add_subplot(gs[1, 3])
cohort_bar(lin_tp, lin_cd,'j gene',ax=jax0,yaxis_upper=0.5, labelsize=12, ticksize=10, legendsize=8)
cohort_bar(seq_tp, seq_cd,'j gene',ax=jax1,yaxis_upper=0.56, labelsize=12, ticksize=10, legendsize=8)
make_djins_plot(lin_tp, jax2, labelsize=12, ticksize=10, legendsize=8)
make_djdel_plot(lin_tp, jax3, labelsize=12, ticksize=10, legendsize=8)
make_djins_plot(seq_tp, jax4, labelsize=12, ticksize=10, legendsize=8)
make_djdel_plot(seq_tp, jax5, labelsize=12, ticksize=10, legendsize=8)
for n, ax in enumerate([jax0, jax1, jax2, jax3, jax4, jax5]):
    ax.text(-0.1, 1.1, "("+string.ascii_uppercase[n]+")", transform=ax.transAxes, 
            size=12, weight='bold')
fig.savefig('/gscratch/stf/zachmon/covid/J_supplementary.pdf',bbox_inches='tight')
