from utils import *
import matplotlib.pyplot as plt

def make_dynamic_bar(time_dict, in_gene, primer=None,
                     ax=None, multiplicity=False):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if in_gene == "v":
        all_genes = abstar_v_genes
    elif in_gene == 'j':
        all_genes = partis_j_genes

    counters = []
    times = list(time_dict.keys())
    for ti in time_dict:
        #if multiplicity == True:
        #    primers = ["IGHV"+str(i) for i in range(1,7)]
        #    for primer in primer:
        #        bad_genes = [gene for gene in time_dict[ti]
        #                     if v_primer not in gene]
        #        for bg in bad_genes:
        #            del time_dict[ti][bg]
        counters.append(get_data_counter(time_dict[ti], multiplicity=multiplicity))

    #  Add all genes
    for counter in counters:
        for g in all_genes:
            if g not in counter:
                counter[g] = 0

    width = 1.0 / len(time_dict) - 0.05
    colors = ['blue', 'orange', 'red', 'black']
    bars = []
    xlabels = []
    default_x = np.arange(0,len(counters[0]),1)
    for i,unsorted_d in enumerate(counters):
        total = sum(unsorted_d.values())
        D = sort_dict(unsorted_d)
        x = default_x+width*i
        xlabels = list(D.keys())
        y = np.divide(list(D.values()),total)
        label_str = "time: " + str(times[i]) + ", # seqs: " + str(total)
        if 'oof' in str(times[i]):
            label_str = str(times[i]) + ", # seqs: " + str(total)
        y = np.divide(list(D.values()),total)
        if i == len(counters) - 1:
            colors[i] = 'black'
        bars.append(ax.bar(x,y,
                           width,
                           color=colors[i]))
                           #label = label_str))


    xtickfontsize = 10
    titlefontsize = 10
    ax.set_yticks(np.arange(0,1+0.1,0.2))

    middle = int(np.mean(np.arange(0,len(counters))))
    ax.set_xticks(default_x + middle*width)
    ax.set_xticklabels(xlabels, fontsize=xtickfontsize, rotation=90)
    #set_loc = 'upper right'
    #if '6' in primer and gene=='v':
    #    set_loc = 'upper left'

def bin_data(data, cdr3=False):
    bin_width = 1
    if cdr3:
        bin_width = 3
    bins = np.arange(0, np.max(data) + bin_width, bin_width)
    data_hist = np.histogram(data, bins)
    data_bins, data_bin_height = data_hist[1][:-1]/bin_width, data_hist[0]
    data_bin_height = np.divide(data_bin_height, len(data))
    return data_bins, data_bin_height


def plot_v_gene_abundance(save_dir, stats):
    gene_type = 'v'
    patient_key = list(stats.keys())[0]
    gene_dict = stats[patient_key]['in frame no indels'][gene_type]
    oof_dict = stats[patient_keys]['out of frame no indels'][gene_type]

    #  Split by primer since using abundance counts
    primer_split = split_dict_by_primer(gene_dict)
    oof_primer_split = split_dict_by_primer(oof_dict)

    #  Split by time
    time_dict = {}
    for i,primer in enumerate(primer_split):
        #  Split ifs by time
        time_dict[primer] = split_dict_by_time(primer_split[primer])
        #  Use oofs from all times as oofs
        time_dict[primer]['oof'] = oof_primer_split[primer]

    fig = plt.figure(dpi=300,figsize=(16,4))
    ax = fig.add_subplot(111)
    ax.margins(0.01)
    ax.grid(True)
    ax.set_axisbelow(True)
    for p in time_dict:
        make_dynamic_bar(patient_key, primer,
                         time_dict[p], gene_type.split(" ")[0],
                         ax=ax, multiplicity=True)

    legend_list = []
    if 'healthy' in patient_key:
        textstr = ("Patient: " + patient_key[index] +
                   "\nSeverity: " + CONST_DATA_DICT[int(patient_key)]['severity'])
        for ti in CONST_DATA_DICT[int(patient_key)]['times']:
            legend_list.append("time="+str(ti))
    else:
        textstr = "Patient: " + patient_key + "\nSeverity: Healthy"
        legend_list.append("time=0")
    legend_list.append("oof")


    ax.legend(legend_list,fontsize=10,loc='upper center')

    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.62, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

    save_name = "patient-"+patient_key+"_" + gene_type.replace(" ", "-") + "_profile_abundance.png"
    fig.savefig(save_dir + save_name, bbox_inches='tight')
    plt.close()

def plot_gene_unique(if_time_dict, oof_dict,gene,
                     ax=None, multiplicity=False,
                     boxtext=None, savename=None):

    #  Split only by time, not primer
    time_dict = if_time_dict
    time_dict['oof'] = oof_dict

    if ax == None:
        fig = plt.figure(dpi=300,figsize=(16,4))
        ax = fig.add_subplot(111)
    ax.margins(0.01)
    ax.grid(True)
    ax.set_axisbelow(True)
    make_dynamic_bar(time_dict, gene,
                     ax=ax, multiplicity=multiplicity)
    legend_list = []
    for ti in time_dict:
        legend_list.append("time="+str(ti))

    ax.legend(legend_list,fontsize=7,loc='upper center')

    if boxtext:
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        if 'j' in gene:
            ax.text(0.1, 0.95, boxtext, transform=ax.transAxes, fontsize=5,
                    verticalalignment='top', bbox=props)
        else:
            ax.text(0.4, 0.95, boxtext, transform=ax.transAxes, fontsize=5,
                    verticalalignment='top', bbox=props)


    if savename:
        fig.savefig(savename, bbox_inches='tight')
        plt.close()

def make_line_scatter_time_dynamic(if_time_dict, oof_dict, data_key, xlabel,
                                   ax=None, multiplicity=False, boxtext=None, savename=None):
    if ax is None:
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)

    if data_key == 'cdr3 length':
        upper_lim = 50
        bins = np.arange(0,3*upper_lim+3,3)
        divisor = 3
    elif 'ins' in data_key:
        upper_lim = 35
        bins = np.arange(0,upper_lim+1,1)
        divisor = 1
    elif 'del' in data_key:
        upper_lim = 25
        bins = np.arange(0,upper_lim+1,1)
        divisor = 1

    if_counters = []
    for ti in if_time_dict:
        if_counters.append(get_data_counter(if_time_dict[ti],multiplicity=multiplicity))

    oof_counters = get_data_counter(oof_dict,multiplicity=multiplicity)

    times = list(if_time_dict.keys())
    max_y_list = []

    #  Oof
    oof_data = []
    for key in oof_counters:
        oof_data += [key]*oof_counters[key]
    num_oof_data = len(oof_data)
    oof_hist = np.histogram(oof_data,bins)
    oof_bins, oof_bin_height = oof_hist[1][:-1]/divisor, oof_hist[0]
    oof_bin_height = np.divide(oof_bin_height, num_oof_data)
    max_y_list.append(np.max(oof_bin_height))
    oof_color = "red"
    ls = ['-']
    ax.plot(oof_bins, oof_bin_height,
            '-', linewidth=1.0, color='black', alpha=1.0,
            marker='X', markersize=2, fillstyle='none',
            label="oof (# seqs=" + str(num_oof_data) + ")")

    colors = ['blue', 'orange', 'red']
    for i, if_counter in enumerate(if_counters):
        ti = str(times[i])
        if_data = []
        for key in if_counter:
            if_data += [key]*if_counter[key]
        if not if_data:
            continue
        num_if_data = len(if_data)

        if_hist = np.histogram(if_data,bins)
        if_bins, if_bin_height = if_hist[1][:-1]/divisor,if_hist[0]
        if_bin_height = np.divide(if_bin_height,num_if_data)
        max_y_list.append(np.max(if_bin_height))
        ls= ['-',':']
        mstyle = ['o','X']
        ax.plot(if_bins, if_bin_height,
                ls[0],linewidth=1.0,color=colors[i],alpha=1.0,
                marker=mstyle[0], markersize=3, fillstyle='none',
                label="ti= " + ti + ", if (# seqs=" + str(num_if_data) + ")")
    ax.set_xlabel(xlabel)
    ax.legend(fontsize=6, loc='upper right')

    max_x = np.max([np.max(oof_bins),np.max(if_bins)])
    max_y = np.max(max_y_list)
    log10_max_y=np.log10(max_y)
    round_y_digit = int(-np.floor(log10_max_y))
    max_y = round(max_y,round_y_digit)
    if round_y_digit+log10_max_y > 0.5:
        y_step = 1*10**(-round_y_digit)
    else:
        y_step = 5*10**(-round_y_digit-1)
    if max_y > 0.7:
        y_step = 0.2
    ax.set_xticks(np.arange(0,max_x+10,10))
    ax.set_yticks(np.arange(0.0,max_y+y_step,y_step))
    ax.set_xlim(0, upper_lim)
    #ax.margins(0.1)

    if boxtext:
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax.text(0.8, 0.4, boxtext, transform=ax.transAxes, fontsize=5,
                verticalalignment='top', bbox=props)

    if savename:
        fig.savefig(savename, bbox_inches='tight')
        plt.close()

def primer_split_plots(if_dict, oof_dict, plotfunc, plotargs,
                       time_threshold=None, multiplicity=False,
                       boxtext=None, savename=None):
    primer_split = split_dict_by_primer(if_dict)
    oof_primer_split = split_dict_by_primer(oof_dict)
    #  Remove entries with multiple primers collapsed
    primer_list = list(primer_split.keys())
    bad_keys = [p for p in primer_list if "," in p]
    for p in bad_keys:
        del primer_split[p]

    num_primers = len(primer_split)
    num_cols = 3
    num_rows = int(num_primers / 3 + 1)
    fig, axs = plt.subplots(num=None, figsize=(16, 8), dpi=300,
                            facecolor='w', edgecolor='k',
                            nrows=num_rows, ncols=num_cols)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=1.0)

    for i,p in enumerate(primer_list):
        if time_threshold:
            time_dict = split_dict_by_time(primer_split[p], time_threshold=time_threshold)['data']
        else:
            time_dict = split_dict_by_time(primer_split[p], time_threshold=time_threshold)
        i_row = int(i/3)
        i_col = int(i%3)
        plotfunc(time_dict, oof_primer_split[p], *plotargs,
                 ax=axs[i_row][i_col],multiplicity=multiplicity,
                 boxtext=boxtext,savename=None)
        axs[i_row][i_col].set_title(p)

    counter = 0
    for ax_row in axs:
        for ax in ax_row:
            counter += 1
            if counter > num_primers:
                fig.delaxes(ax)

    if savename:
        fig.savefig(savename, bbox_inches='tight')
        plt.close()

def make_pvalue_scatter_plots(fisher_output, ax=None):
    if ax is None:
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)

    colors = ["blue", "orange", "red", "green", "purple", "deeppink"]
    markers = ["o","s","p","X",'D',"*"]
    fisher_output = sort_dict(fisher_output)
    for j, primer in enumerate(fisher_output):
        num_clusters = len(fisher_output[primer]['early'])
        # Get number of points that will show on plot
        time_diverse_clusters = sum((fisher_output[primer]['early'] != 0)
                                    * (fisher_output[primer]['late'] != 0))
        legend_label = ("primer " + primer +
                       ", total # clusters=" + str(num_clusters) +
                       "(#shown="+str(time_diverse_clusters)+")")
        ax.scatter(np.log10(fisher_output[primer]['oddsratio']),
                   np.log10(fisher_output[primer]['pvalue']),
                   marker=markers[j],color=colors[j],
                   alpha=1.0,s=7)
        ax.scatter(np.log10(fisher_output[primer]['oddsratio'][fisher_output[primer]['pvalue'] == 0.0]),
                   [-340]*len(fisher_output[primer]['oddsratio'][fisher_output[primer]['pvalue'] == 0.0]),
                   marker=markers[j],color=colors[j],
                   alpha=1.0,s=7,label=legend_label)
    ax.legend(fontsize=5, loc='lower left')
    ax.set_xlabel("log10[oddsratio]")
    ax.set_ylabel("log10[pvalues]")
    ax.set_xlim(-5,5)
    ax.set_ylim(-350,5)
    ax.legend(fontsize=5)
    return ax

def make_expansion_plots(fisher_output, count_type, ax=None):
    if ax is None:
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)

    colors = ["blue", "orange", "red", "green", "purple", "deeppink"]
    markers = ["o","s","p","X",'D',"*"]
    fisher_output = sort_dict(fisher_output)
    for j, primer in enumerate(fisher_output):
        num_clusters = len(fisher_output[primer]['early'])
        normalized_early = (fisher_output[primer]['early']
                            / (fisher_output[primer]['early'] + fisher_output[primer]['other early']))
        normalized_late = (fisher_output[primer]['late']
                           / (fisher_output[primer]['late'] + fisher_output[primer]['other late']))

        # Get number of points that will show on plot
        time_diverse_clusters = sum((fisher_output[primer]['early'] != 0)
                                    * (fisher_output[primer]['late'] != 0))

        legend_label = ("primer " + primer +
                        ", total # clusters=" + str(num_clusters) +
                        "(#shown="+str(time_diverse_clusters)+")")
        pvalue_thresh = 1e-200
        ax.scatter(np.log10(normalized_early[fisher_output[primer]['pvalue'] < pvalue_thresh]),
                   np.log10(normalized_late[fisher_output[primer]['pvalue'] < pvalue_thresh]),
                    marker=markers[j], color=colors[j],
                    alpha=1.0, s=7, label=legend_label)
        ax.scatter(np.log10(normalized_early[fisher_output[primer]['pvalue'] > pvalue_thresh]),
                   np.log10(normalized_late[fisher_output[primer]['pvalue'] > pvalue_thresh]),
                    marker=markers[j], color=colors[j],
                    alpha=0.5, s=7, facecolors='None', linewidth=0.4)
    ax.legend(fontsize=5, loc='upper left')
    ax.set_xlabel("log10[early " + count_type + " fraction]")
    ax.set_ylabel("log10[late " + count_type + " fraction]")
    ax.set_xlim(-6,0)
    ax.set_ylim(-6,0)
    ax.plot(ax.get_xlim(),
            ax.get_ylim(),
            ls="--", c=".3")
    ax.legend(fontsize=5)
    return ax
