#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Tiago Oliveira and D. Pratas"
__copyright__ = ("Copyright (C) 2020, IEETA, University of Aveiro."
                 "This is a Free software, under GPLv3. You may redistribute"
                 "copies of it under the terms of the GNU - General Public"
                 "License v3 <http://www.gnu.org/licenses/gpl.html>. There"
                 "is NOT ANY WARRANTY, to the extent permitted by law.")
__credits__ = ["Tiago Oliveira","Diogo Pratas"]
__license__ = "GPL-3.0"
__version__ = "DLCR - 0.1"
__maintainer__ = "Tiago Oliveira"
__email__ = "tiagomanuel28@gmail.com"
__status__ = "Production"
__description__ = "DLCR: Efficient detection of Distant Low Complexity Regions "


## Imports

import csv
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import tqdm
from matplotlib import gridspec, ticker

import src.file_cleanup as cleanup


def check_fasta(fasta):
    '''
    Checks if fasta file path is correct
    :param fasta: fasta file path
    :return: correct fasta path
    '''
    fasta = str(fasta)
    check = cleanup.check_folder_change("src", "root", fasta)
    if check:
        return str(check)
    else:
        return fasta


def check_positions(fasta):
    '''
    Checks positions file path is correct
    :param fasta: fasta file path
    :return: correct positions file path
    '''
    positions = str(fasta) + ".positions"
    check = cleanup.check_folder_change("src", "root", fasta)
    if check:
        return str(check)
    else:
        return positions


def get_poisitions(fasta):
    '''
    Extracts sequences presentes in positions file from fasta file, using gto_fasta_extract
    :param fasta: fasta file path
    :return: list of tupples, containing seq and starting position
    '''
    sequences = []
    positions = check_positions(fasta)
    fasta = check_fasta(fasta)
    with open(positions, "r") as positions_file:
        with open(fasta, "r") as fasta_file:
            for line in positions_file:
                if len(line) > 0:
                    seq = ""
                    begin = int(line.split(":")[0])
                    end = int(line.split(":")[1])
                    command = ["gto_fasta_extract", "-i", str(begin), "-e", str(end)]
                    seq = subprocess.run(args=command, stdin=fasta_file, stdout=subprocess.PIPE, check=True)
                    seq = seq.stdout.decode('utf-8')
                    if len(seq) >= 20:
                        sequences.append((seq,begin))
    return sequences

def get_reps(fasta):
    '''
    Gets low complexity regions from .positions file, calculates de length and returns a list with the values
    :param fasta: name of the FASTA file profiled for complex regions
    :return: list with length values
    '''
    rep_sizes = []
    positions = check_positions(fasta)
    with open(positions, "r") as positions_file:
        for line in positions_file:
            if len(line) > 0:
                size = int(line.split(":")[1]) - int(line.split(":")[0])
                if size >= 20:
                    rep_sizes.append(size)

    return rep_sizes

def get_profile(file):
    '''
    Gets the pair of coordinates and compression value of profiling low complexity regions
    :return: tuple with lists, with coordinates and compression value
    '''
    file = check_fasta(file)
    coordinates = list()
    compression = list()
    with open(file, "r") as compression_result:
        for line in compression_result:
            if len(line) > 0:
                coordinates.append(int(line.split("\t")[0]))
                compression.append(float(line.split("\t")[1]))

    return coordinates, compression

def profile_draw_graphs(ax1, ax2, rep_sizes, hist_max, name, threshold, nested=False, bar_space=None, minor_size=4):
    '''
    Draws subplots for profile reporting, given 2 subplots (ax1 and ax2)
    :param ax1: Subplot from a GridSpec Object. Used to draw compression profile plot
    :param ax2: Subplot from a GridSpec Object. Used to draw histogram with the length frequency of low
    complexity regions
    :param rep_sizes: list with lengths low complexity regions
    :param hist_max: int with maximum value of y axis to be used in ax2
    :param name: str with name of the file to be saved in .pdf format
    :param threshold: float with threshold used during profiling, to be drawn as a horizontal line ax1
    :param nested: Bol with information if the GridSpec object is nested
    :param bar_space: float with space to be used for separation of bars in histogram (ax2)
    :param minor_size: int with font size of minor label in x axis of histogram (ax2)
    '''
    # Get compression profile information
    coordinates, compression = get_profile("PROFILE_N")
    x_scale = 10 ** (len(str(len(coordinates))) - 1)
    x_limit = (int(len(coordinates) / x_scale) + 1) * x_scale

    # Compression profile (plot)
    ax1.plot(coordinates, compression, linewidth=0.5, color="b")
    ax1.plot([0, x_limit], [threshold, threshold], 'k-', lw=0.5, label="_not in legend", color="g")
    ax1.set_xlim([0,x_limit])
    ax1.set_ylim([0, 2])
    ax1.grid(True, alpha=0.5)
    ax1.set_ylabel("Bps")
    ax1.set_xlabel=("Length (Bps)")
    ax1.set_title("Low Complexity Regions")

    # Histogram
    if len(rep_sizes) == 0:
        ax2.hist(rep_sizes, bins=1, align="mid", rwidth=bar_space, edgecolor='k')
    else:
        ax2.hist(rep_sizes, bins=len(rep_sizes), align="mid", rwidth=bar_space, edgecolor='k')
    ax2.set_ylabel("Frequency of Regions Length")
    ax2.set_xlabel("Length (Bps)")
    ax2.set_title("Low Complexity Regions Sizes")
    ax2.grid(axis='y', alpha=0.5)
    ax2.set_ylim([0,hist_max+1])
    ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.tick_params(which='minor', length=1, labelsize=minor_size)
    ax2.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%d'))

    # Save figure in .pdf format
    if not nested:
        plt.savefig(name)



def profile_fasta(fasta, level, window, threshold, drop):
    '''
    Runs modified gto_profile_regions.sh for rfms
    :param fasta: str with name of the fasta file to profiled
    :param level: int with level to be used by gto in profiling
    :param window: int with length of the window to be used in profiling
    :param threshold: float with threshold to be used to acess compression (value between 0 and 2)
    :param drop: int with drop (sensibility of compression)
    :return:
    '''
    #TODO: NEED TO ADD TRY METHOD IN CASE PROFILING FAILS DURING .SH RUN
    command = ["sh", "src/gto_rfms_profile_regions.sh", fasta, str(level), str(window), str(threshold), str(drop)]
    subprocess.call(command, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

def profile_report(fasta, level, window, threshold, drop, multiple=False, show=False, positions=True):
    '''
    Functions generates a report for the profiling of low complexity regions in a genomic sequence using gto
    compression tools
    :param fasta: str with name of the .FASTA file to profiled
    :param level: int with level to be used internally in GTO for compression profiling
    :param window: int with length of window to be used for compression profiling
    :param threshold: float with threshold for compression profiling (between 0.0 and 2.0)
    :param drop: int with drop of compression profiling (sensibility of compression)
    :param multiple: bol with information if a multiple window report should be done (default=False)
    :param show: bol with information if the graph should be displayed in a window after completion (default=False)
    :param positions bol indicates if positions file should be cleaned after report completion
    '''
    cleanup.create_folder("graphs")
    fasta_name = str(fasta.split(".")[0] +
                     "_l" + str(level) +
                     "_w" + str(window) +
                     "_t" + str(threshold) +
                     "_d" + str(drop))

    profile_graph = str("graphs/" + fasta_name) + ".pdf"
    profile_report_name = str("graphs/" + fasta_name) + "_report.pdf"

    if multiple: # Multiple window report
        # Setting up tqdm progress bar
        bar = tqdm.tqdm(total=5)
        # Initialization of figure with size of an A4 page
        f = plt.figure(figsize=(8.27, 11.69))
        plt.rc('axes', labelsize=6)
        plt.rc('axes', titlesize=6)  # fontsize of the axes title
        plt.rc('xtick', labelsize=4)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=4)  # fontsize of the tick labels
        # Creation of a main GridSpec object
        main_grid = gridspec.GridSpec(15,7, wspace=0.4, hspace=1.5)

        ## Report for window / 4
        # Nesting of GridSpec object
        grid1 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[:4,:3])
        #Handling of progress bars
        description = str("Generating report for Window: " + str(int(window) / 4))
        bar.set_description_str(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(int(window)/4), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax1 = plt.subplot(grid1[:2, 0:])
        ax2 = plt.subplot(grid1[2:, 0:])
        ax1.text(0.5, 1.3, str("Window: " + str(int(window)/4)), size=10, ha="center", transform=ax1.transAxes)
        #Drawing of subplots in figure
        profile_draw_graphs(ax1, ax2, rep_sizes, hist.max(), profile_report_name, float(threshold), True, 0.5, 2)
        # Progress update
        bar.update()

        ## Report for window / 2
        # Nesting of GridSpec object
        grid2 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[:4, 4:])
        # Handling of progress bars
        description = str("Generating report for Window: "+ str(int(window) / 2))
        bar.set_description(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(int(window) / 2), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax3 = plt.subplot(grid2[:2, 0:])
        ax4 = plt.subplot(grid2[2:, 0:])
        ax3.text(0.5, 1.3, str("Window: " + str(int(window)/2)), size=10, ha="center", transform=ax3.transAxes)
        # Drawing of subplots in figure
        profile_draw_graphs(ax3, ax4, rep_sizes, hist.max(), profile_report_name, float(threshold), True, 0.5, 2)
        # Progress update
        bar.update()

        ## Report for window
        # Nesting of GridSpec object
        grid3 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[5:9, 2:5])
        # Handling of progress bars
        description = str("Generating report for Window: "+ str(int(window)))
        bar.set_description_str(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(window), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax5 = plt.subplot(grid3[:2, 0:])
        ax6 = plt.subplot(grid3[2:, 0:])
        ax5.text(0.5, 1.3, str("Window: " + str(int(window))), size=10, ha="center", transform=ax5.transAxes)
        # Drawing of subplots in figure
        profile_draw_graphs(ax5, ax6, rep_sizes, hist.max(), profile_report_name, float(threshold), True, 0.5, 2)
        # Progress update
        bar.update()

        ## Report for window * 2
        # Nesting of GridSpec object
        grid4 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[10:, :3])
        # Handling of progress bars
        description = str("Generating report for Window: " + str(int(window) * 2))
        bar.set_description_str(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(int(window) * 2), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax7 = plt.subplot(grid4[:2, 0:])
        ax8 = plt.subplot(grid4[2:, 0:])
        ax7.text(0.5, 1.3, str("Window: " + str(int(window)*2)), size=10, ha="center", transform=ax7.transAxes)
        # Drawing of subplots in figure
        profile_draw_graphs(ax7, ax8, rep_sizes, hist.max(), profile_report_name, float(threshold), True, 0.5, 2)
        # Progress update
        bar.update()

        ## Report for window * 4
        # Nesting of GridSpec object
        grid5 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[10:, 4:])
        # Handling of progress bars
        description = str("Generating report for Window: " + str(int(window) * 4))
        bar.set_description(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(int(window) * 4), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax9 = plt.subplot(grid5[:2, 0:])
        ax10 = plt.subplot(grid5[2:, 0:])
        ax9.text(0.5, 1.3, str("Window: " + str(int(window)*4)), size=10, ha="center", transform=ax9.transAxes)
        # Drawing of subplots in figure
        profile_draw_graphs(ax9, ax10, rep_sizes, hist.max(), profile_report_name, float(threshold), False, 0.5, 2)
        # Handling of progress bars
        bar.update()
        description = "Finished compiling profile report"
        bar.set_description(desc=description)

    else: #If not nested
        # Setting up tqdm progress bar
        bar = tqdm.tqdm(total=2, desc="Generating complexity profile")
        # Profiling
        profile_fasta(fasta, level, window, threshold, drop) # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Handling of progress bars
        bar.set_description_str(desc="Generating profile report")
        bar.update()
        # Creation of GridSpec object
        grid = plt.GridSpec(4, 3, wspace=0.4, hspace=2)
        #title="Complexity profile of " + str(fasta) + \
        #      "\nlevel: " + str(level) + \
        #      ", window: " + str(window) + \
        #      ", threshold" + str(threshold) + \
        #      ", drop: " + str(drop)
        #plt.suptitle(title, fontsize=12)
        # Creation of Subplots
        ax1 = plt.subplot(grid[:2, 0:])
        ax2 = plt.subplot(grid[2:, 0:])
        # Drawing of subplots in figure
        profile_draw_graphs(ax1, ax2, rep_sizes, hist.max(), profile_graph, float(threshold), multiple, show)
        # Handling of progress bars
        bar.set_description_str(desc="Finished compiling profile report")
        bar.update()

    cleanup.clean_profiling(positions) # Cleanup of files produced during profiling
    if show: # Show figure in a seperate window
        plt.show()

def reverse_seq(seq):
    '''
    Reverses DNA sequence
    :param seq: string of DNA sequence
    :return: string with reverse DNA sequence
    '''
    bases = {"A":"T",
             "T":"A",
             "C":"G",
             "G":"C"}
    new_seq=""
    for base in seq:
        new_seq += bases[base.upper()]
    return new_seq

def regions_info(fasta):
    '''
    Generates dictionary with information about the low complexity regions
    it contains the sequence, initial position, reversed sequence probability, related regions,
    minimum tandem pieces and maximum tandem pieces
    :param fasta: fasta file path
    :return: dict with low complexity region information
    '''
    seqs = {}
    i_seqs = 0
    for seq in get_poisitions(fasta):
        seqs[i_seqs] = {"seq":seq[0], "pos":seq[1], "reversed":0, "related":[], "min_tandem":0, "max_tandem":0}
        i_seqs += 1
    return seqs

def seq_compression(seq1, seq2=""):
    '''
    Compresses 2 sequences. If only 1 sequence is giver it compresses only 1
    :param seq1: string with seq1 to compress
    :param seq2: string with seq2 to compress
    '''
    if len(seq2) > len(seq1):
        temp_seq = seq1
        seq1 = seq2
        seq2 = temp_seq
    seq = str(seq1).upper()+str(seq2).upper()

    with open("comparison.seq", "w") as comparison_seq:
        comparison_seq.write(seq)
    command = ["src/DLCR", "-F", "comparison.seq"]
    #seq = subprocess.run(args=command, stdout=subprocess.PIPE, check=True) #prints running time
    seq = subprocess.run(args=command, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL, check=True)

def move_search_cursor(array, pointer):
    '''
    moves search cursor in tandem search when switching between down threshold and up threshold
    :param array: list of indexes of values relative to threshold
    :param pointer: int with cursor
    :return: tuple containing shortned list of indexes and updated search cursor
    '''
    array = array[array > pointer]
    if array.size > 0:
        return array, int(array[0])
    else:
        return array, int(pointer+1)

def tandem_compression_profile(drop = 20, threshold = 1.75, compression = 0.70):
    '''
    Compresses each sequence and determines the possible ammount of tandems present
    in sequence.
    :param drop: sensibility of search in bps
    :param threshold: threshold of the compression
    :param compression: threshold to determine compression
    :return: ints with information about size and number of tandems present in sequence
    '''
    reps_size = []
    profile1 = np.loadtxt(".dlcr_1.inf", delimiter=";")
    profile2 = np.loadtxt(".dlcr_2.inf", delimiter=";")
    profile2 = profile2[::-1] # .dlcr_2.dna in reverse order for both side profiling
    comparison = profile1 <= profile2
    profile = profile2
    for i in range(0, profile1.size):
        if comparison[i]:
            profile[i] = profile1[i]

    if profile.size < 60:
        drop = int(profile.size/3)

    compression_profile = profile[profile <= float(threshold)].size
    compression_profile = compression_profile / profile.size
    if compression_profile > float(compression):
        profile3 = np.where(profile2 <= threshold)[0]
        i = profile3[0]
        for j in np.where(profile2 > threshold)[0]:
            if j-i >= drop and (j-i) >= 4:
                reps_size.append(j-i)
        i = profile3[profile3 > j]
        profile1 = np.where(profile <= threshold)[0]
        profile2 = np.where(profile > threshold)[0]
    else:
        profile1 = np.where(profile <= threshold)[0]
        profile2 = np.where(profile > threshold)[0]
        if profile1.size == 0:
            return False, False, False, False
        else:
            i = profile1[0]
    down_thresh = True
    start = False
    while i + drop < profile.size:
        sliced_profile = profile[i:(i+drop+1)]
        compression_profile = sliced_profile[sliced_profile <=float(threshold)].size
        if compression_profile/sliced_profile.size > float(compression):
            if down_thresh:
                if not start:
                    start = i
                profile2, i = move_search_cursor(profile2, i)
            else: # down_thresh
                if (i-start) <= 4:
                    profile1, i = move_search_cursor(profile1, i)
                else:
                    reps_size.append(i-start)
                    profile1, i = move_search_cursor(profile1, i)
                    start = False

        else: # compression threshold
            if down_thresh:
                if start:
                    profile2, i = move_search_cursor(profile2, i)
                else: #start
                    profile2, i = move_search_cursor(profile2, i)
                    profile1, i = move_search_cursor(profile1, i)
                    down_thresh = False
            else: # down_thresh
                profile1, i = move_search_cursor(profile1, i)
        # preparing new cycle
        if down_thresh:
            down_thresh = False
        else:
            down_thresh = True
    #outside while
    if start:
        end = profile.size
        reps_size.append(end - start)
    else:
        reps_size.append(profile.size)
    size_min = min(reps_size)
    size_max = max(reps_size)
    rep_min = int(profile.size/size_max)
    rep_max = int(profile.size/size_min)
    return size_min, size_max, rep_min, rep_max

def region_compression_profile(seq1_size, threshold = 1.75, compression = 0.75):
    '''
    Acesses if seq2 is compressed with seq1
    :param seq1_size: size of seq1
    :param threshold: log2 complexity score for relation threshold
    :param compression: relation compression rassio
    :return: Bolean indicating if seq2 is present in seq1
    '''
    profile = np.loadtxt(".dlcr_1.inf", delimiter=";")[seq1_size:]
    if (profile[profile <= float(threshold)].size / profile.size) > float(compression):
        return True # Same
    else:
        return False # Different

def regions_associate_known(low_complexity_seqs, threshold=1.75, compression = 0.70):
    #load known sequences from file
    known_seqs = {}
    for key in low_complexity_seqs.keys():
        for known_seq in known_seqs.keys():
            seq_compression(known_seqs[known_seq] ,low_complexity_seqs[key]["seq"])
            if region_compression_profile(len(known_seqs[known_seq]), threshold,
                                          compression):
                pass

def seq_reference(low_complexity_seqs, key):
    '''
    Function generates seq reference to use in FALCON
    :param low_complexity_seqs: dic with low complexity regions
    :param key: key of low_complexity_seq to be used as ref
    '''
    with open("ref_seq.fasta", "w") as ref_seq:
        fasta_header = ">" + str(key) + "\n"
        ref_seq.write(fasta_header)
        ref_seq.write(low_complexity_seqs[key]["seq"])

def seq_database(low_complexity_seqs, keys):
    '''
    Function generates seq_database to use in FALCON
    :param low_complexity_seqs: dic with low complexity regions
    :param keys: list of low_complexity_seq keys that are related
    '''
    with open("seq_db.fasta", "w") as seq_db:
        for key in keys:
            fasta_header = ">"+str(key)+"\n"
            seq_db.write(fasta_header)
            seq_db.write(low_complexity_seqs[key]["seq"])
            seq_db.write("\n\n")

def similarity_read(relation_compression):
    '''
    Function that parses top.txt file, resulting of FALCON processing
    :param relation_compression: compression % used to determine relationship between regions
    :return: list of similar seqs that correspond to low_complexity_seqs keys
    '''
    similar = []
    with open("top.txt") as top:
        for line in top:
            compression =  float(line.replace("\n", "").split("\t")[2])
            if compression >= relation_compression:
                similar.append(int(line.replace("\n", "").split("\t")[3]))
            else:
                break
    return similar

def run_falcon(level=str(4), top_size=str(20), threads=str(2)):
    '''
    Run Falcon for similaraity assessment
    :param level: int with FALCON level
    :param top_size: number of top seqs for top.txt
    :param threads: int of threads for precessin
    '''
    threads = str(threads)
    command = ["FALCON", "-F", "-l", level, "-t", top_size, "-n", threads, "-x", "top.txt",
               "ref_seq.fasta", "seq_db.fasta"]
    seq = subprocess.run(args=command, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, check=True)

def regions_assotiate_reversed(low_complexity_seqs, key, related, relation_compression, threads=str(2)):
    '''
    function determines if low complexity regions are associated by reverse sequence
    :param low_complexity_seqs: dic with low complexity regions
    :param key: key of low_complexity_seq to be used as ref
    :param related: list of low_complexity_seq keys that are related
    :param relation_compression: compression % used to determine relationship between regions
    :param threads:  int of threads for precessing
    :return: dic with low complexity regions
    '''

    seq_database(low_complexity_seqs, related)
    level = str(3)
    top_size = str(len(related))
    threads = str(threads)
    run_falcon(level, top_size, threads)
    similar = similarity_read(relation_compression)
    if similar:
        low_complexity_seqs[key]["reversed"] = float(len(similar)/len(related))
    return low_complexity_seqs


def regions_associate_different(low_complexity_seqs, relation_compression, threads=2):
    '''
    function determines relation between low complexity regions
    :param low_complexity_seqs: dic with low complexity regions
    :param relation_compression: compression % used to determine relationship between regions
    :param threads: int of threads for precessing
    :return: dic with low complexity regions
    '''
    regions_keys = list(low_complexity_seqs.keys())
    level = str(4)
    top_size = str(len(regions_keys))
    threads = str(threads)

    with tqdm.tqdm(total=len(regions_keys), desc="relations", position=2, leave=False) as pbar:
        for i_key in range(0,len(regions_keys)):
            seq_reference(low_complexity_seqs, regions_keys[i_key])
            seq_database(low_complexity_seqs, regions_keys[i_key+1:])
            run_falcon(level, top_size, threads)
            similar = similarity_read(relation_compression)
            if similar:
                temp_similar = low_complexity_seqs[regions_keys[i_key]]["related"]
                temp_similar = temp_similar + similar
                low_complexity_seqs[regions_keys[i_key]]["related"] = temp_similar
                for similar_region in similar:
                    temp_similar = low_complexity_seqs[similar_region]["related"]
                    temp_similar += [regions_keys[i_key]]
                    low_complexity_seqs[similar_region]["related"] = temp_similar
                low_complexity_seqs = regions_assotiate_reversed(low_complexity_seqs, regions_keys[i_key],
                                                                 similar, relation_compression, threads)
            pbar.update(1)
    return low_complexity_seqs


def regions_associate_tandems(low_complexity_seqs, drop=20, tandem_threshold=1.75, tandem_compression=0.7):
    '''
    Determines tandem sequences present inside low complexity regions identified low complexity regions in report
    :param low_complexity_seqs: dict with low complexity region information
    :param drop: drop to be used during tandem seq search
    :param tandem_threshold: log2 complexity threshold for tandem compression
    :param tandem_compression: compression % used to determine presence of tandem sequence
    :return: dict with low complexity regions information
    '''
    with tqdm.tqdm(total=len(list(low_complexity_seqs.keys())), desc="Managing Tandems", position=2, leave=False) as tbar:
        for key in low_complexity_seqs.keys():
            seq_compression(low_complexity_seqs[key]["seq"])  # Compress single sequence
            size_min, size_max, rep_min, rep_max = tandem_compression_profile(drop, tandem_threshold, tandem_compression)
            if size_min:
                low_complexity_seqs[key]["min_tandem"] = rep_min
                low_complexity_seqs[key]["max_tandem"] = rep_max
            tbar.update(1)
    return low_complexity_seqs

def regions_associate(low_complexity_seqs, drop = 20, tandem_threshold=1.75, tandem_compression=0.7,
                      relation_compression=0.75, threads=2):
    '''
    Determines association between different sequences, tandem and low complexity regions
    :param low_complexity_seqs: dict with low complexity region information
    :param drop: drop to be used during tandem seq search
    :param tandem_threshold: log2 complexity threshold for tandem compression
    :param tandem_compression: compression % used to determine presence of tandem sequence
    :param relation_compression: compression % used to determine relation between 2 regions
    :param threads: int with number of threads for relation profile processing
    :return: dict with low complexity regions information
    '''
    #TODO:ADD known sequences
    #low_complexity_seqs = regions_associate_known(low_complexity_seqs, tandem_threshold, tandem_compression)
    low_complexity_seqs = regions_associate_tandems(low_complexity_seqs, drop, tandem_threshold, tandem_compression)

    keys = list(low_complexity_seqs.keys())
    low_complexity_seqs = regions_associate_different(low_complexity_seqs, relation_compression, threads=2)

    return low_complexity_seqs


def calculate_percentages(fasta):
    '''
    Calculates sequence information to create profile INI for RFMS
    :param fasta: fasta file path
    :return: returns ACTG % and seq size or False if problem when reading Fasta File
    '''
    fasta = check_fasta(fasta)
    with open(fasta, "r") as file:
        passed = False
        seq = ""
        a = 0
        t = 0
        g = 0
        c = 0
        size = 0
        for line in file:
            if not line.startswith("<") and line.replace("\n",""):
                seq = str(line).replace("\n", "").lower()
                a += seq.count("a")
                t += seq.count("t")
                g += seq.count("g")
                c += seq.count("c")
                size += len(seq)
        if size > 0:
            a = float(a / size)
            t = float(t / size)
            g = float(g / size)
            c = float(c / size)
            return a, t, g, c, size
        else:
            return False, False, False, False, False



def generate_INI(fasta, output_file, org_name, sample_type):
    '''
    Generates profile INI file to be used in RFMS
    :param fasta: fasta file path
    :param output_file: output file path
    :param org_name: org name to be atributed to INI
    :param sample_type: sample type information
    :return: size of full sequence
    '''
    a, t, g, c, size = calculate_percentages(fasta)
    if a:
        with open(output_file, "a") as file:
            file.write("["+str(org_name)+"]\n")
            file.write("a="+str(a)+"\n")
            file.write("t="+str(t)+"\n")
            file.write("g="+str(g)+"\n")
            file.write("c="+str(c)+"\n")
            file.write("seq_min="+str(size)+"\n")
            file.write("seq_max="+str(size+1)+"\n")
            file.write("sample_type="+str(sample_type)+"\n")
            return size
    else:
        raise OSError(str(fasta) + " is in invalid FASTA format")


def generate_csv_dict(low_complexity_seqs, size_seq):
    '''
    Generates dict used to create CSV file with low complexity region information
    :param low_complexity_seqs: dict with low complexity regions information
    :param size_seq: size full sequence
    :return: dict with information for csv writing
    '''
    csv_dict = {}
    regions = list(low_complexity_seqs.keys())
    regions_parsed = set()
    repeat = 1
    for region in regions:
        related = low_complexity_seqs[region]["related"].copy()
        if region not in regions_parsed or (region in regions_parsed and related):
            nr_regions = 1
            size = np.array(len(low_complexity_seqs[region]["seq"]))
            min_tandem = np.array(int(low_complexity_seqs[region]["min_tandem"]))
            max_tandem = np.array(int(low_complexity_seqs[region]["max_tandem"]))
            reversed = np.array(float(low_complexity_seqs[region]["reversed"]))
            windows = {}
            pos = int(low_complexity_seqs[region]["pos"])
            window_size = int(size_seq/10)
            for window in  range(1,11):
                if pos < (window_size*window) and pos >= ((window_size*window)-window_size):
                    windows[window] = [pos]
                else:
                    windows[window] = []
            for related_region in related:
                nr_regions+=1
                min_tandem = np.append(min_tandem, int(low_complexity_seqs[related_region]["min_tandem"]))
                max_tandem = np.append(max_tandem, int(low_complexity_seqs[related_region]["max_tandem"]))
                reversed = np.append(reversed, float(low_complexity_seqs[related_region]["reversed"]))
                low_complexity_seqs[region]["related"].remove(related_region)
                low_complexity_seqs[related_region]["related"].remove(region)
                size = np.append(size, np.array(len(low_complexity_seqs[related_region]["seq"])))
                for i in related:
                    if related_region in low_complexity_seqs[i]["related"]:
                        low_complexity_seqs[i]["related"].remove(related_region)
                pos = int(low_complexity_seqs[related_region]["pos"])
                for window in range(1, 11):
                    if pos < (window_size * window) and pos >= ((window_size * window) - window_size):
                        windows[window].append(pos)

            regions_parsed.update(related)
            min_tandem = float(min_tandem.min())
            max_tandem = float(max_tandem.max())
            min_size = int(size.min())
            max_size = int(size.max())
            if min_size == max_size:
                max_size += 1
            if nr_regions == 1:
                min_nr_regions = 1
                max_nr_regions = 2
            else:
                min_nr_regions = int(nr_regions - (nr_regions / 2))
                max_nr_regions = int(nr_regions + (nr_regions / 2))
            reversed = round(float(reversed.mean()),3)
            window1_prob = str(round(len(windows[1])/nr_regions, 3))
            if windows[1]:
                window1_std = str(round(np.array(windows[1]).std(),3))
                window1_mean = str(round(np.array(windows[1]).mean(),3))
                if window1_std == "0.0":
                    window1_std = str(round(window_size*0.05,3))
            else:
                window1_std = str(0.0)
                window1_mean = str(0.0)
            window2_prob = str(round(len(windows[2])/nr_regions,3))
            if windows[2]:
                window2_std = str(round(np.array(windows[2]).std(),3))
                window2_mean = str(round(np.array(windows[2]).mean(),3))
                if window2_std == "0.0":
                    window2_std = str(round(window_size*0.05,3))
            else:
                window2_std = str(0.0)
                window2_mean = str(0.0)
            window3_prob = str(round(len(windows[3])/nr_regions,3))
            if windows[3]:
                window3_std = str(round(np.array(windows[3]).std(),3))
                window3_mean = str(round(np.array(windows[3]).mean(),3))
                if window3_std == "0.0":
                    window3_std = str(round(window_size*0.05,3))
            else:
                window3_std = str(0.0)
                window3_mean = str(0.0)
            window4_prob = str(round(len(windows[4])/nr_regions,3))
            if windows[4]:
                window4_std = str(round(np.array(windows[4]).std(),3))
                window4_mean = str(round(np.array(windows[4]).mean(),3))
                if window4_std == "0.0":
                    window4_std = str(round(window_size*0.05,3))
            else:
                window4_std = str(0.0)
                window4_mean = str(0.0)
            window5_prob = str(round(len(windows[5])/nr_regions,3))
            if windows[5]:
                window5_std = str(round(np.array(windows[5]).std(),3))
                window5_mean = str(round(np.array(windows[5]).mean(),3))
                if window5_std == "0.0":
                    window5_std = str(round(window_size*0.05,3))
            else:
                window5_std = str(0.0)
                window5_mean = str(0.0)
            window6_prob = str(round(len(windows[6])/nr_regions,3))
            if windows[6]:
                window6_std = str(round(np.array(windows[6]).std(),3))
                window6_mean = str(round(np.array(windows[6]).mean(),3))
                if window6_std == "0.0":
                    window6_std = str(round(window_size*0.05,3))
            else:
                window6_std = str(0.0)
                window6_mean = str(0.0)
            window7_prob = str(round(len(windows[7])/nr_regions,3))
            if windows[7]:
                window7_std = str(round(np.array(windows[7]).std(),3))
                window7_mean = str(round(np.array(windows[7]).mean(),3))
                if window7_std == "0.0":
                    window7_std = str(round(window_size*0.05,3))
            else:
                window7_std = str(0.0)
                window7_mean = str(0.0)
            window8_prob = str(round(len(windows[8])/nr_regions,3))
            if windows[8]:
                window8_std = str(round(np.array(windows[8]).std(),3))
                window8_mean = str(round(np.array(windows[8]).mean(),3))
                if window8_std == "0.0":
                    window8_std = str(round(window_size*0.05,3))
            else:
                window8_std = str(0.0)
                window8_mean = str(0.0)
            window9_prob = str(round(len(windows[9])/nr_regions,3))
            if windows[9]:
                window9_std = str(round(np.array(windows[9]).std(),3))
                window9_mean = str(round(np.array(windows[9]).mean(),3))
                if window9_std == "0.0":
                    window9_std = str(round(window_size*0.05,3))
            else:
                window9_std = str(0.0)
                window9_mean = str(0.0)
            window10_prob = str(round(len(windows[10])/nr_regions,3))
            if windows[10]:
                window10_std = str(round(np.array(windows[10]).std(),3))
                window10_mean = str(round(np.array(windows[10]).mean(),3))
                if window10_std == "0.0":
                    window10_std = str(round(window_size*0.05,3))
            else:
                window10_std = str(0.0)
                window10_mean = str(0.0)


            csv_dict[repeat] = {"min_tandem":min_tandem, "max_tandem":max_tandem, "min_size":min_size,"max_size":max_size,
                                "min_nr_regions":min_nr_regions, "max_nr_regions":max_nr_regions, "reversed":reversed,
                                "window1":str(window1_prob+":"+window1_mean+":"+window1_std),
                                "window2":str(window2_prob+":"+window2_mean+":"+window2_std),
                                "window3":str(window3_prob+":"+window3_mean+":"+window3_std),
                                "window4":str(window4_prob+":"+window4_mean+":"+window4_std),
                                "window5":str(window5_prob+":"+window5_mean+":"+window5_std),
                                "window6":str(window6_prob+":"+window6_mean+":"+window6_std),
                                "window7":str(window7_prob+":"+window7_mean+":"+window7_std),
                                "window8":str(window8_prob+":"+window8_mean+":"+window8_std),
                                "window9":str(window9_prob+":"+window9_mean+":"+window9_std),
                                "window10":str(window10_prob+":"+window10_mean+":"+window10_std),
                                "repeat":repeat}
            repeat+=1
    return csv_dict


def generate_csv(org_name, csv_dict):
    '''
    Writes csv file with information present in csv_dict, to be used in RFMS
    :param org_name: org name to be used in filename
    :param csv_dict: dict with low complexity regions information ready to write
    '''
    filename = str(org_name) + ".csv"
    with open(filename, 'w') as csvfile:
        fieldnames = ["repeat", "min_tandem", "max_tandem", "min_size", "max_size", "min_nr_regions", "max_nr_regions",
                      "reversed", "window1", "window2", "window3", "window4", "window5", "window6", "window7",
                      "window8", "window9", "window10"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=";")
        writer.writeheader()
        for key in list(csv_dict.keys()):
            writer.writerow(csv_dict[key])


def generate_profile(fasta, level, window, threshold, drop, output, org_name, sample_type,
                     drop_tandem = 20, tt=1.75, tc=0.7, rc=0.75, threads=2, show=False):
    '''
    Functions generates a profile of low complexity regions in a genomic sequence using gto
    compression tools and FALCON.
    :param fasta: str with name of the .FASTA file to profiled
    :param level: int with level to be used internally in GTO for compression profiling
    :param window: int with length of window to be used for compression profiling
    :param threshold: float with threshold for compression profiling (between 0.0 and 2.0)
    :param drop: int with drop of compression profiling (sensibility of compression)
    :param output: str of output INI file
    :param org_name: str with org name for INI file and profile
    :param sample_type: str with sample type of org
    :param drop_tandem: int with drop to be used during tandem seq search
    :param tt: float log2 complexity threshold for tandem compression
    :param tc: int compression % used to determine presence of tandem sequence
    :param rc: int compression % used to determine relation between 2 regions
    :param threads: int number of threads for relation profiling processing
    :param show: bol with information if the graph should be displayed in a window after completion (default=False)
    :return:
    '''

    with tqdm.tqdm(total=6, desc="Generating Profiles", position=1, leave=True) as mbar:
        if show:
            profile_report(fasta, level, window, threshold, drop, False, show, False)
        else:
            profile_fasta(fasta, level, window, threshold, drop)
            cleanup.clean_profiling(False)
        mbar.update(1)


        low_complexity_seqs = regions_info(fasta)
        mbar.update(1)
        low_complexity_seqs = regions_associate(low_complexity_seqs, drop_tandem, tt, tc, rc, threads)
        cleanup.clean_profiling(True)
        cleanup.clean_dlcr()
        mbar.update(1)
        size = generate_INI(fasta, output, org_name, sample_type)
        mbar.update(1)
        #generate csv
        mbar.update(1)
        csv_dict = generate_csv_dict(low_complexity_seqs, size)
        generate_csv(org_name, csv_dict)
        mbar.update(1)


# Testing and debugging
if __name__ == '__main__':
    pass
    #profile_report("BS_genome_NC_000964.3.fasta", 5, 10042156, 1.75, 0, False, True)
    generate_profile("BS_genome_NC_000964.3.fasta", 5, 10042156, 1.75, 0, False, False)