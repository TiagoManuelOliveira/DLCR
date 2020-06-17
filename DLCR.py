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
import argparse
import configparser  # use to store profiles in INI format
from datetime import datetime

from src.generate_profiles import generate_profile
from src.generate_profiles import profile_report


def greetings():
    greetings_text = "                                    \n" \
                     "                                    \n" \
                     "######   ##        #####   #####    \n" \
                     "#######  ##       #######  #######  \n" \
                     "##   ##  ##       ##            ##  \n" \
                     "##   ##  ##       ##       ######   \n" \
                     "##   ##  ##       ##       ## ##    \n" \
                     "#######  #######  #######  ##  ##   \n" \
                     "######    ######   #####   ##   ##  \n" \
                     "                                    \n" \
                     "                                    "

    print(greetings_text)

def report_func(args):
    print("report")

def profile_func(args):
    print("profile")

#TODO: Add epilogues with usage for argparse?
#TODO: Add copywrite in help?

# Generate argparser objects
parser = argparse.ArgumentParser(prog="DLCR",description=__version__ + " - " + __description__, epilog=greetings())
parser.add_argument("-v", "--version", action="version", version=__version__)
# Generate subparsers for report and profile operations
subparser = parser.add_subparsers(help="sub-command help", required=True, dest="mode")

report = subparser.add_parser("report", help="Creates reports of distant low complexity regions",
                              description="Creates reports of distant low complexity regions")
profile = subparser.add_parser("profile", help="Create profile of distant low complexity regions for generation "
                                               "of simulated reads in RFMS",
                               description="Create profile of distant low complexity regions for generation of "
                                           "simulated reads in RFMS")
# Arguments for report generation
report.set_defaults(func=report_func)
report_files = report.add_argument_group(title="file arguments")
report_variables=report.add_argument_group(title="variable arguments")
report_display = report.add_argument_group(title="report display arguments")
report_files.add_argument("-f", "--fasta", action="store", help="Input file name for FASTA file", default=False,
                        required=True, metavar="file.fasta")

report_variables.add_argument("-l", "--level", action="store", help="Level for profiling", default=False,
                        required=True, metavar=3)
report_variables.add_argument("-w", "--window", action="store", help="window length for complexity profile", default=1000,
                    required=False, metavar=1000)
report_variables.add_argument("-t", "--threshold", action="store", help="complexity log2 threshold", default=1.5,
                    required=False, metavar=1.5)
report_variables.add_argument("-d", "--drop", action="store", help="drop of bps in complexity profile", default=0,
                    required=False, metavar=0)
report_display.add_argument("-m", "--multiple", action="store_true", help="multiple window report", required=False)
report_display.add_argument("-s", "--show", action="store_true", help="display report", required=False)

# Arguments for profile generation
profile.set_defaults(func=profile_func)
profile_files = profile.add_argument_group(title="file arguments")
profiles_variables = profile.add_argument_group(title="variable arguments")
profile_display = profile.add_argument_group(title="display arguments")
profile_files.add_argument("-f", "--fasta", action="store", help="Input file name for FASTA file", default=False,
                        required=True, metavar="file.fasta")
profiles_variables.add_argument("-l", "--level", action="store", help="Level for profiling", default=False,
                        required=True, metavar=3)
profiles_variables.add_argument("-w", "--window", action="store", help="window length for complexity profile",
                                default=1000, required=False, metavar=1000)
profiles_variables.add_argument("-t", "--threshold", action="store", help="complexity log2 threshold", default=1.5,
                     required=False, metavar=1.5)
profiles_variables.add_argument("-d", "--drop", action="store", help="drop of bps in complexity profile", default=0,
                     required=False, metavar=0)
profiles_variables.add_argument("-tt", "--tandem_threshold", action="store", help="complexity log2 threshold for tandem profiling",
                     default=1.75, required=False, metavar=1.75, dest="tt")
profiles_variables.add_argument("-tc", "--tandem_compression", action="store", help="compression threshold for tandem profiling",
                     default=0.7, required=False, metavar=0.7, dest="tc")
profiles_variables.add_argument("-td", "--tandem_drop", action="store", help="Drop used during tandem profiling",
                                default=20, required=False, metavar=20, dest="drop_tandem")
profiles_variables.add_argument("-rt", "--relation_threads", action="store", help="relation compression number of threads",
                     default=2, required=False, metavar=2, dest="threads")
profiles_variables.add_argument("-rc", "--relation_compression", action="store", dest="rc",
                                help="compression threshold for relation profiling", default=0.75, required=False, metavar=0.75)
profile_display.add_argument("-s", "--show", action="store_true", help="generates, saves and displays report of "
                                                                       "profile", required=False)
profile_files.add_argument("-o", "--output", action="store", help="Output file name for profile", required=False,
                     metavar="<output>.profile", default=False, dest="output")
profile_files.add_argument("-n", "--org_name", help="Output organism name", action="store", required=False,
                           metavar="ORG", dest="org_name", default="ORG")
profile_files.add_argument("-st", "--sample_type", help="type of organism", action="store", required=False,
                           metavar="bacteria", dest="sample_type", default="NA")
#print(parser.argument_default)
args = parser.parse_args()

## Files
profiles_file = "INIs/profiles.ini"
orgs_file = "INIs/orgs.ini"
profiles = configparser.ConfigParser()
profiles.read(profiles_file)
orgs = configparser.ConfigParser()
orgs.read(orgs_file)

## main

def arg_common(args):
    '''
    Parse arguments common to both modes
    :param args: argparse object
    :return: returns common arguments
    '''
    # fasta
    fasta = str(args.fasta)
    # level
    try:
        level = int(args.level)
    except:
        raise ValueError("Level must be an integer")
    # window
    #TODO:Windowsize = len/1000 como default
    try:
        window = int(args.window)
    except:
        raise ValueError("Window must be an integer")
    # threshold
    try:
        threshold = float(args.threshold)
        if threshold > 2 or threshold < 0:
            raise ValueError("threshold value must be between 2 and 0")
    except:
        raise ValueError("threshold must be a float")
    # drop
    try:
        drop = int(args.drop)
        if drop < 0:
            raise ValueError("Drop must be greater than 0")
        if drop > window:
            raise ValueError("Drop can't be greater than window")
    except:
        raise ValueError("drop must be an integer")
    # show
    show = args.show


    return fasta, level , window, threshold, drop, show

def arg_profile(args):
    '''
    Parse arguments specific to profile mode
    :param args: argparse object
    :return: returns parsed profile arguments and returns report arguments as False
    '''
    # multiple
    multiple = False
    # output
    if args.output:
        output = str(str(args.output) + ".profile")
    else:
        current_date = str(datetime.date(datetime.now())).replace("-", "_")
        current_time = "_".join(str(datetime.time(datetime.now())).split(":")[:2])
        output = str("DLCR_Profile_" + current_time + "_" + current_date + ".profile")
    # org_name
    org_name = str(args.org_name)
    # tandem_threshold
    try:
        tt = float(args.tt)
        if tt > 2 or tt < 0:
            raise ValueError("tandem threshold value must be between 2 and 0")
    except:
        raise ValueError("tandem threshold must be a float")
    # tandem compression
    try:
        tc = float(args.tc)
        if tc > 1 or tc < 0:
            raise ValueError("tandem compression value must be between 1 and 0")
    except:
        raise ValueError("tandem compression must be a float")
    # relation compression
    try:
        rc = float(args.rc)
        if rc > 1 or rc < 0:
            raise ValueError("relation compression value must be between 1 and 0")
    except:
        raise ValueError("relation compression must be a float")
    # Threads
    try:
        threads = int(args.threads)
    except:
        raise ValueError("Thread must be integer")
    # sample type
    sample_type = str(args.sample_type)
    # tandem drop
    drop_tandem = int(args.drop_tandem)

    return multiple, output, org_name, tt, tc, rc, threads, sample_type, drop_tandem

def arg_report(args):
    '''
    Parse arguments specific to report mode
    :param args: argparse object
    :return: returns specific report arguments and returns profile arguments as False
    '''
    # multiple
    multiple = args.multiple
    # output
    output = False
    # org_name
    org_name = False
    # tandem_threshold
    tt = False
    # tandem compression
    tc = False
    # relation compression
    rc = False
    # Threads
    threads = False
    # sample type
    sample_type = False
    # tandem drop
    drop_tandem = False

    return multiple, output, org_name, tt, tc, rc, threads, sample_type, drop_tandem


def arg_handler():
    '''
    Handles arguments from argparse
    :return: returns all alguments from arparse
    '''
    #mode
    mode = str(args.mode)

    fasta, level, window, threshold, drop, show = arg_common(args)

    if mode == "profile":
        multiple, output, org_name, tt, tc, rc, threads, sample_type, drop_tandem = arg_profile(args)

    if mode == "report":
        multiple, output, org_name, tt, tc, rc, threads, sample_type, drop_tandem = arg_report(args)

    return mode, fasta, level, window, threshold, drop, multiple, show, output, org_name, tt, tc, rc, threads, sample_type, drop_tandem

def run():
    '''
    Runs DLCR:
    Calls function to parse arguments in argparse and then runs one of 2 modes: report or profile
    report mode, creates low complexity reports in pdf format
    profile mode, analyzes low complexity regions and produces a profile compliant with RFMS
    :return:
    '''
    mode, fasta, level, window, threshold, drop, multiple, show, output, org_name, tt, tc, rc, threads, sample_type, drop_tandem = arg_handler()
    if mode == "report":
        #TODO Modify this prints in final version
        print(mode)
        print(fasta, level, window, threshold, drop, multiple, show)
        profile_report(fasta, level, window, threshold, drop, multiple, show)
    else:
        print(mode)
        print(fasta, level, window, threshold, drop, output, org_name, tt, tc, rc, threads, show)
        generate_profile(fasta, level, window, threshold, drop, output, org_name, sample_type, drop_tandem, tt, tc, rc, threads, show)

if __name__ == '__main__':
    run()