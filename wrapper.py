#!/usr/bin/env python3

import subprocess
import glob
import shutil
import os


# Function that concatenates all the .fastq.gz files in each individual barcode folder and renames them after
# their respective folder

def concatenate(rootdir):
    print("Concatenating files")    
    folder_names = glob.glob(f"{rootdir}/*")
   
    for i in range(0, len(folder_names)):   # loop through all barcode folders
        folder = os.path.basename(os.path.normpath(folder_names[i]))    # Extract the name of the folder (e.g. barcode01)
        bc_file = glob.glob(f"{rootdir}{folder}/*.fastq.gz")    # barcode file with full path
        with open(f"{folder}.fastq.gz",'wb') as wfd:
            for f in bc_file:
                with open(f,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)    # this copies to pwd (where the script is run from). How do I copy to a specific folder?

concatenate("path/") # Input path to the folder containing the barcode folders