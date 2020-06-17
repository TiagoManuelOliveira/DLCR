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

import os

def create_folder(folder):
    '''
    If folder for project doesnt exist, create it
    :param folder: folder name
    '''
    if os.path.exists(folder) == False:
        print("Folder ", folder, " doesnt exist...\n>", folder, " folder created")
        os.makedirs(folder)

def clean_profiling(positions=True):
    '''
    Cleans files created during reporting/general profiling
    :param positions:Bol, keep positions file
    '''
    files = ["A_D",
             "A_min",
             "A_R",
             "FIL_UB.x",
             "FIL_UB_N.x",
             "FIL_UB_R.x",
             "IDXES",
             "PROFILE_N",
             "SEQ",
             "SEQ.co",
             "SEQ.iae",
             "SEQ_R",
             "SEQ_R.co",
             "SEQ_r.iae",
             "SEQ_R.iae",
             "SEQ_R_UB",
             "SEQ_UB"]

    for file in os.listdir():
        # find .positions files
        if file[-10:] == ".positions":
            if positions:
                files.append(file)
        # find .seq files
        if file[-4:] == ".seq":
            files.append(file)


    for file in files:
        if file in os.listdir():
            os.remove(file)

def clean_dlcr():
    '''
    Clean temporary profiling files
    '''
    files = [".dlcr_1.dna",
             ".dlcr_2.dna",
             ".dlcr_1.inf",
             ".dlcr_2.inf",
             "comparison.seq",
             "ref_seq.fasta",
             "seq_db.fasta",
             "top.txt",
             "local.fal"]

    for file in files:
        if file in os.listdir():
            os.remove(file)


def check_folder_change(origin, target, file):
    '''
    Check if file exists in origin folder. If not changes file path to target
    :param origin: str origin folder path
    :param target: str target folder path
    :param file: str file path to change
    :return: False if file in origin folder path. New path if file in origin folder path
    '''
    try:
        file = str(file)
    except:
        raise TypeError("filename must be type str")
    if file in os.listdir():
        return False
    else:
        if origin == "src" and target == "root":
            if file in os.listdir(".."):
                return "../" + file
            else:
                raise OSError(file +" can't be found")


if __name__ == '__main__':
    clean_profiling()