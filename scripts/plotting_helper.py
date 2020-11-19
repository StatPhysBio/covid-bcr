#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script containing helper functions for creating plots and performing statistical tests.
    Copyright (C) 2020 Montague, Zachary

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from math import log10, floor
import string

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import numpy as np
from scipy.stats import f_oneway, ks_2samp
import seaborn as sns

from utils import *

#  dictionary used when referencing times, severity, plotting binding, etc.
#CONST_DATA_DICT = get_bcell_info('/gscratch/stf/zachmon/covid/plasma_data/plasma_info.csv')#total_bcell/old_total_b_cell_info.csv')
CONST_DATA_DICT = get_bcell_info('/gscratch/stf/zachmon/covid/total_bcell/old_total_b_cell_info.csv')
PLASMA_DATA_DICT = get_bcell_info('/gscratch/stf/zachmon/covid/plasma_data/plasma_info.csv')
for key in PLASMA_DATA_DICT:
    CONST_DATA_DICT[key] = PLASMA_DATA_DICT[key]
#get_bcell_info('/gscratch/stf/zachmon/covid/plasma_data/plasma_info.csv')

#  Set fonttypes so that Adobe Illustrator can be used to edit the final product.
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

#  Set colorscheme for data.
colors = {'Healthy':"#0072B2",
          'Mild': "#009E73",
          'Moderate':"#E69F00",
          "Severe": "#D55E00",
          'Asymptomatic': "#CC79A7"}

#  Set lighter colors to distinguish from darker colors.
lightercolors = {'Healthy': (144/255, 215/255, 255/255),
                 'Mild': (46/255, 255/255, 200/255),
                 'Moderate': (255/255, 212/255, 110/255),
                 "Severe": (255/255, 138/255, 47/255),
                 'Asymptomatic': (255/255, 26/255, 255/255)}

#  Function to round to n significant digits.
round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))

def geo_mean(iterable: np.array) -> np.array:
    """Returns the geometric mean of an iterable.

    Parameters
    ----------
    iterable : np.array
        Array containing values.

    Returns
    -------
    np.float64
        Geometric mean of the array.
    """

    return np.exp(np.mean(np.log(iterable)))

def geo_std(iterable: np.array) -> np.array:
    """Returns the geometric standard deviation of an iterable.

    Parameters
    ----------
    iterable : np.array
        Array containing values.

    Returns
    -------
    np.float64
        Geometric mean of the array.
    """

    return np.exp(np.std(np.log(iterable)))

def bin_and_average_values(counters: list, observable: str, indict: bool = False,
                           geo: bool = False) -> (np.array, np.array, list):
    """Bins the values within a cohort and obtains the average and variation within that bin.

    Parameters
    ----------
    counters : list
        List of dictionaries containing the count information for each individual
        in a cohort. The keys are effectively the fine-grained bins whereas the values
        are the counts in each fine-grained bin.
    observable : str
        The quantity which is to be binned and averaged.
    indict : bool, optional
        Specifies if the input is a dictionary (for genes) or list (all other observables).
    geo :  bool, optional
        If True, will perform geometric averaging, otherwise arithmetic averaging.

    Returns
    -------
    avgs : np.array
        Array of the averages in each bin.
    stds : np.array
        Array of one standard deviation of variation in each bin.
    x : list
        List of the values of bins which might be strings (for genes) or integers.
    """

    #  Binning for observables.
    binning = {'cdr3 length': np.linspace(2,36,18),
                 'vd ins': np.linspace(0,30,31),
                 'dj ins': np.linspace(0,30,31),
                 'vd del': np.linspace(0,35,36),
                 'dj del': np.linspace(0,42,43)}

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

def cohort_plot(cohort_averages: dict, observable: str, in_ax: matplotlib.axes = None,
                geo: bool = False) -> matplotlib.axes:
    """Creates a plot with a line for the average values of a cohort and shading for the variation.

    Parameters
    ----------
    cohort_averages : dict
        Dictionary of averages and variations within in a cohort for all severities.
    observable : str
        The quantity which is going to be plotted.
    in_ax : matplotlib.axes, optional
        The axes on which values are going to be plotted.
        Used to modify an already existing axes when creating a figure with a grid.
    geo : bool, optional
        Specifies whether the variation should be +/- or *//.

    Returns
    -------
    ax : matplotlib.axes
        Axes on which the plot is plotted.
    """

    if in_ax is None:
        fig, ax = plt.subplots(figsize=(4, 4),dpi=300)
    else:
        ax = in_ax
    for s in cohort_averages:
        x = cohort_averages[s][observable][-1]
        ax.plot(x,
                cohort_averages[s][observable][0],
                alpha=1.0,
                color=colors[s], label=s, linewidth = 1.0)
        if geo:
            ax.fill_between(x,
                            cohort_averages[s][observable][0]/cohort_averages[s][observable][1],
                            cohort_averages[s][observable][0]*cohort_averages[s][observable][1],
                            color=colors[s], alpha=0.15, linewidth=0)
        else:
            ax.fill_between(x,
                            cohort_averages[s][observable][0]-cohort_averages[s][observable][1],
                            cohort_averages[s][observable][0]+cohort_averages[s][observable][1],
                            color=colors[s], alpha=0.15, linewidth=0)
    if in_ax is None:
        return fig,ax
    else:
        return ax

def format_axes(in_ax: matplotlib.axes, labelsize: int = 12, ticksize: int = 10,
               legendsize: int = 8, box: bool = False) -> None:
    """Formats ticks, tick label sizes, and axes label sizes for box plots.

    Parameters
    ----------
    in_ax : matplotlib.axes
        Axes containing a plot.
    labelsize : int, optional
        Size of axes labels.
    ticksize : int, optional
        Size of tick labels.
    legendsize : int, optional
        Specifies font size of strings in legend.
    box : bool, optional
        Specifies whether or not a box plot has been plotted.

    Returns
    -------
    None
    """
    in_ax.xaxis.set_minor_locator(AutoMinorLocator())
    in_ax.yaxis.set_minor_locator(AutoMinorLocator())
    in_ax.xaxis.label.set_size(labelsize)
    in_ax.yaxis.label.set_size(labelsize)
    in_ax.tick_params(axis='both', which='major', direction='in',length=4,
                      bottom=True, top=True, left=True, right=True,
                      labelsize=ticksize, color='grey')
    in_ax.tick_params(axis='both', which='minor', direction='in',length=2,
                      bottom=True, top=True, left=True, right=True,
                      color='grey')
    if not box:
        in_ax.legend(fontsize=legendsize)

def make_stats_plot(cohort_averages: dict, observable: str, in_ax: matplotlib.axes,
                    geo: bool = False, labelsize: int = 12, ticksize: int = 10,
                    legendsize: int = 8, box=False) -> None:
    """Plots an observable for all cohorts on a single axes object.

    Parameters
    ----------
    cohort_averages : dict
        Dictionary of averages and variations within a cohort for all severities.
    observable : str
        The quantity which is going to be plotted.
    in_ax : matplotlib.axes
        The axes on which values are going to be plotted.
    geo: bool, optional
        Specifies geometric or arithmetic averaging.
    labelsize : int, optional
        Size of axes labels.
    ticksize : int, optional
        Size of tick labels.
    legendsize : int, optional
        Specifies font size of strings in legend.
    box : bool, optional
        Specifies whether or not a box plot has been plotted.

    Returns
    -------
    None
    """

    cohort_plot(cohort_averages, observable, in_ax=in_ax, geo=geo)

    #  Name of xlabel for each observable.
    xlabels = {'cdr3 length': 'HCDR3 length [aa]',
               'vd ins': 'VD insertions [nt]',
               'vd del': 'VD deletions [nt]',
               'dj ins': 'DJ insertions [nt]',
               'dj del': 'DJ deletions [nt]'}
    in_ax.set_xlabel(xlabels[observable],fontname="Arial")
    in_ax.set_ylabel("PDF", family="Arial")
    format_axes(in_ax, labelsize=labelsize, ticksize=ticksize,
                legendsize=legendsize, box=False)

def cohort_bar(cohort_averages: dict, cohort_data: dict, observable: str,
               yaxis_upper: int = None, ax: matplotlib.axes = None, labelsize: int =12,
               ticksize: int = 10, legendsize: int = 8) -> None:
    """Plots either the V- or J-gene usages.

    Parameters
    ----------
    cohort_averages : dict
        Dictionary of averages and variations within a cohort for all severities.
    cohort_dict : dict
        Dictionary of all statistics of all individuals in a cohort for all severities.
    observable : str
        The quantity which is going to be plotted.
    yaxis_upper : int, optional
        Specifies the upper limit of the y-axis on the plot.
    ax : matplotlib.axes, optional
        Used to modify an already existing axes when creating a figure with a grid.
    labelsize : int, optional
        Size of axes labels.
    ticksize : int, optional
        Size of tick labels.
    legendsize : int, optional
        Specifies font size of strings in legend.

    Returns
    -------
    None
    """

    if ax is None:
        fig = plt.figure(dpi=300,figsize=(16,4))
        ax = fig.add_subplot(111)

    #  Give the bars for each gene some space among each other.
    width = 1.0 / len(cohort_averages) - 0.02
    bars = []

    #  Sort the genes by descending usage in the healthy cohort.
    sorted_genes = [gene
                    for _,gene in sorted(zip(cohort_averages['Healthy'][observable][0],
                                             cohort_averages['Healthy'][observable][-1]),
                                              key=lambda pair: pair[0], reverse=True)]

    #  Because there are so many V genes, plot only those which have at least
    #  1% average usage in at least one cohort.
    if 'v' in observable:
        good_genes = []
        for gene in sorted_genes:
            bools = 0
            for severity in cohort_averages:
                idx = cohort_averages[severity][observable][-1].index(gene)
                value = cohort_averages[severity][observable][0][idx]
                bools += value >= 0.01
            if bools != 0:
                good_genes.append(gene)
    else:
        good_genes = sorted_genes

    xlabels = good_genes
    default_x = np.arange(0,len(good_genes),1)
    xs,ys = [], []

    for i,severity in enumerate(cohort_averages):
        x = default_x + width*i
        xs.append(x)
        indices = [cohort_averages[severity][observable][-1].index(gene)
                   for gene in good_genes]
        y = [cohort_averages[severity][observable][0][k] for k in indices]
        ys.append(y)

    if observable == 'v gene':
        ax.set_xlim(-0.5, xs[0][-1] + 1.0)
        if yaxis_upper is not None:
            ax.set_ylim(0, yaxis_upper)
        else:
            ax.set_ylim(0, 0.35)
    else:
        ax.set_xlim(-0.3, xs[0][-1] + 1.0)
        if yaxis_upper is not None:
            ax.set_ylim(0, yaxis_upper)
        else:
            ax.set_ylim(0, 0.55)

    #  Specifiy location of xticks as being in the middle of the grouping of bars.
    middle = np.mean(np.arange(0, len(cohort_averages)))
    middle_xticks = default_x + middle*width
    ax.set_xticks(middle_xticks)
    ax.set_xticklabels(xlabels, rotation=90, family='Arial')
    ax.set_ylabel('PDF', fontname="Arial")

    #  Plot the bars. 
    for i,severity in enumerate(cohort_averages):
        bar = ax.bar(xs[i], ys[i], width, color=colors[severity], label=severity)
        for d in cohort_data[severity][observable]:
            #  Plot the full distributions of values to get a better sense
            #  of variation within cohorts.
            for gidx,rect in enumerate(bar):
                ax.scatter(xs[i][gidx], d[good_genes[gidx]], marker='.',
                           color=lightercolors[severity], zorder=3, s=15,
                           edgecolors='black', linewidths=0.3)

    format_axes(ax, labelsize=labelsize, ticksize=ticksize, legendsize=legendsize)

def get_stats(in_files, rep_type='bulk', statstype='progenitors', geo=False, pseudocount=0) -> (dict, dict):
    """Obtains the statistics from each patient and puts them into their respective cohort and also obtains averaged cohort statistics.

    Parameters
    ----------
    in_files : list
        List of json files containing statistics output.
    statstype: str, optional
        Options: progenitors, nonsingletons.
    geo : bool, optional
        Specifies whether geometric averaging or arithmetic averaging is performed.
    pseudocount: int, optional
        A pseudocount of 0 is used for arithmetic averaging whereas a value of 1
        is used in geometric averaging. Be very careful using pseudocounts when
        performing statistics tests, particularly ANOVA.

    Returns
    -------
    cohort_data : dict
        Dictionary of all data for all individuals in a cohort for all severities.
    cohort_averages : dict
        Dictionary of average data within in cohort for all severities.
    """

    stats = {}
    for f in in_files:
        instats = json_open(f)[statstype]
        patient = instats['patient'][0]
        stats[patient] = {'v gene':{}, 'j gene':{}, 'cdr3 length': {},
                                          'vd ins':{}, 'dj ins': {},
                                          'vd del':{}, 'dj del': {}}
        for observable in stats[patient]:
            try:
                temp_list = [o for idx, o in enumerate(instats[observable])
                             if instats[rep_type][idx]]
                stats[patient][observable] = Counter(temp_list)
            except:
                stats[patient][observable] = Counter(instats[observable])
        #  Leave out empty stats
        if not list(stats[patient][observable].values()):
            del stats[patient]

    #  Names of observables to access from statistics output.
    observables = ['v gene', 'j gene', 'cdr3 length', 'vd ins',
                   'dj ins', 'vd del', 'dj del']

    cohort_raw_data = {}
    cohort_data = {}
    cohort_averages = {}
    for key in colors:
        cohort_raw_data[key] = []
        cohort_data[key] = {}
        cohort_averages[key] = {}

    for patient in stats:
        severity = CONST_DATA_DICT[patient]['severity']
        if str(severity) == 'nan':
            severity = 'Mild'
        for o in observables:
            cohort_data[severity][o] = []
            cohort_averages[severity][o] = []
        cohort_raw_data[severity].append(stats[patient])

    for severity in colors:
        if not cohort_data[severity]:
            del cohort_data[severity]
            del cohort_averages[severity]
            del cohort_raw_data[severity]

    for severity in cohort_raw_data:
        for observable in cohort_raw_data[severity][0]:
            counters = [ind_stats[observable]
                        for ind_stats in cohort_raw_data[severity]]
            equalize_counters(counters, pseudocount)
            for counter in counters:
                if observable == 'v gene' or observable == 'j gene':
                    cohort_data[severity][observable].append(sort_dict_by_value(counter))
                else:
                    datalist = []
                    for value in counter:
                        datalist += [value]*counter[value]
                    #  Convert CDR3 length from nt length to aa length
                    if observable == 'cdr3 length':
                        datalist = np.array(datalist) / 3
                    cohort_data[severity][observable].append(datalist)
            if observable == 'v gene' or observable == 'j gene':
                continue
            else:
                cohort_averages[severity][observable] = bin_and_average_values(cohort_data[severity][observable],
                                                                            observable, indict=False, geo=False)

    equalize_counters([counter for severity in cohort_data
                       for counter in cohort_data[severity]['v gene']], pseudocount)
    equalize_counters([counter for severity in cohort_data
                       for counter in cohort_data[severity]['j gene']], pseudocount)

    for severity in cohort_data:
        for observable in ['v gene', 'j gene']:
            for counter in cohort_data[severity][observable]:
                normalize_counter(counter)
    for severity in cohort_averages:
        for observable in ['v gene', 'j gene']:
            cohort_averages[severity][observable] = bin_and_average_values(cohort_data[severity][observable],
                                                                       observable, indict=True, geo=False)
    return cohort_data, cohort_averages

def perform_anova(cohort_data: dict, observable: str, ax=None) -> None:
    """Creates box plots and performs ANOVA on the means of each individual's statistic.

    Parameters
    ----------
    cohort_data : dict
        Dictionary of all data for all individuals in a cohort for all severities.
    observable : str
        The quantity which is going to be plotted and tested.

    Returns
    -------
    None
    """

    #  Name of ylabel for each observable.
    ylabels = {'cdr3 length': 'HCDR3 length [aa]',
               'vd ins': 'VD insertions [nt]',
               'vd del': 'VD deletions [nt]',
               'dj ins': 'DJ insertions [nt]',
               'dj del': 'DJ deletions [nt]'}

    means = {}
    for severity in cohort_data:
        means[severity] = []
        for ind_data in cohort_data[severity][observable]:
            means[severity].append(np.mean(ind_data))

    for severity in means:
        if severity == 'Healthy':
            continue
        f_val, p_val = f_oneway(means['Healthy'], means[severity])
        df_groups = 1
        df_sample_size = len(means['Healthy']) + len(means[severity]) - 2
        print('Healthy vs. ' + severity + ':',
              'F_{'+ str(df_groups) + ',' + str(df_sample_size) + '} =',
              str(round_to_n(f_val, 3)) + ',',
              'p-value =',
              round_to_n(p_val, 3))

    if ax is None:
        fig=plt.figure(dpi=300)
        ax = fig.add_subplot(111)

    sns.set(context='paper', style='white')
    plt.rcParams.update({"xtick.bottom" : True, "ytick.left" : True})
    ax = sns.boxplot(data=list(means.values()), width=.9, palette=colors.values(), ax=ax)
    sns.swarmplot(data=list(means.values()), size=6, edgecolor="black",
                  linewidth=.9, palette=lightercolors.values(), ax=ax)
    ax.set_xlabel('Cohort',fontname='Arial')
    ax.set_xticklabels(list(means.keys()), fontname='Arial')
    ax.set_ylabel('Mean ' + ylabels[observable], fontname='Arial')
    yticklabels = ax.get_yticks()
    ax.set_yticklabels(yticklabels, fontname='Arial')
    format_axes(ax, box=True)
    plt.show()

def perform_ks2samp(cohort_data: dict, observable: str) -> None:
    """Performs two-sample K–S test on pooled cohort data.

    Parameters
    ----------
    cohort_data : dict
        Dictionary of all data for all individuals in a cohort for all severities.
    observable : str
        The quantity which is going to be tested.

    Returns
    -------
    None
    """

    for severity in cohort_data:
        if severity == 'Healthy':
            continue
        all_healthy = np.concatenate([ind_data
                                      for ind_data in cohort_data['Healthy'][observable]])
        single_cohort = np.concatenate([ind_data
                                        for ind_data in cohort_data[severity][observable]])
        ksstat,pval = ks_2samp(all_healthy, single_cohort)
        healthy_len = len(all_healthy)
        cohort_len = len(single_cohort)
        print('Healthy vs. ' + severity + ':',
              'D_{' + str(healthy_len) + ',' + str(cohort_len) +'} =',
              str(round_to_n(ksstat,3)) + ',',
              'p value =',round_to_n(pval,3))

